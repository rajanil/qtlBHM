import random
import cPickle
import gzip
import sys
import pdb
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import preprocess

def plot_annotation_weights(weights, stderrs, labels):

    figure = plot.figure(figsize=(7,10))
    subplot = figure.add_axes([0.1,0.1,0.8,0.8])
    subplot.spines['top'].set_color('none')
    subplot.spines['right'].set_color('none')
    subplot.spines['left'].set_color('none')
    subplot.yaxis.set_ticks([])
    subplot.xaxis.set_ticks_position('bottom')

    # select significant annotations
    sig = np.logical_or((weights-stderrs)>0,(weights+stderrs)<0)
    weights = weights[sig]
    stderrs = stderrs[sig]
    labels = [labels[i] for i,s in enumerate(sig) if s]

    # order annotations
    order = np.argsort(weights)
    weights = weights[order]
    stderrs = stderrs[order]
    labels = [labels[o] for o in order]

    N = order.size
    colormask = np.zeros((N,), dtype='int')
    enriched = (weights-stderrs)>0
    depleted = (weights+stderrs)<0
    colormask[depleted] = -1
    colormask[enriched] = 1

    yvals = np.arange(1,2*N,2)
    xmin = np.min(weights-stderrs)-0.5
    xmax = np.max(weights+stderrs)

    # significantly enriched
    if np.any(enriched):
        subplot.errorbar(weights[enriched], yvals[enriched], xerr=stderrs[enriched], marker='o',\
            capsize=0, linestyle='none', color='r', elinewidth=1, markeredgewidth=1)

    # significantly depleted
    if np.any(depleted):
        subplot.errorbar(weights[depleted], yvals[depleted], xerr=stderrs[depleted], marker='o', \
            capsize=0, linestyle='none', color='b', elinewidth=1, markeredgewidth=1)

    for i,m in enumerate(colormask):
        if m==-1:
            plot.text(xmin, yvals[i], labels[i], fontsize=8, color='b', \
                horizontalalignment='right', verticalalignment='center')
        elif m==1:
            plot.text(xmin, yvals[i], labels[i], fontsize=8, color='r', \
                horizontalalignment='right', verticalalignment='center')

    subplot.axis([xmin-5, xmax, -1, 2*N+1])
    subplot.axvline(0, linestyle='--', color='k', linewidth=1)
    
    subplot.axes.get_yaxis().set_visible(False)

    subplot.set_xlabel('annotation weights', fontsize=10)
    for tick in subplot.get_xticklabels():
        tick.set_fontsize(10)

    return figure

def plot_annotation_proportion(preprop, postprop, labels):

    figure = plot.figure(figsize=(7,10))
    subplot = figure.add_axes([0.1,0.1,0.8,0.8])
    subplot.spines['top'].set_color('none')
    subplot.spines['right'].set_color('none')
    subplot.spines['left'].set_color('none')
    subplot.yaxis.set_ticks([])
    subplot.xaxis.set_ticks_position('bottom')

    order = np.argsort(postprop-preprop)
    if order.size>20:
        order = np.hstack((order[:10],order[-10:]))
    preprop = preprop[order]
    postprop = postprop[order]
    labels = [labels[o] for o in order]

    N = order.size
    yvals = np.arange(1,2*N,2)
    xmin = 0
    xmax = max([np.max(postprop), np.max(preprop)])+0.1

    enriched = postprop>preprop
    depleted = postprop<=preprop

    # significantly enriched
    if np.any(enriched):
        subplot.scatter(postprop[enriched], yvals[enriched], marker='o', \
            color='r', linewidth=1)
        subplot.scatter(preprop[enriched], yvals[enriched], marker='o', \
            color='#888888', linewidth=1)

    # significantly depleted
    if np.any(depleted):
        subplot.scatter(postprop[depleted], yvals[depleted], marker='o', \
            color='b', linewidth=1)
        subplot.scatter(preprop[depleted], yvals[depleted], marker='o', \
            color='#888888', linewidth=1)

    for i,m in enumerate(enriched):
        if m:
            plot.text(xmin-0.25, yvals[i], labels[i], fontsize=8, color='r', \
                horizontalalignment='right', verticalalignment='center')
        else:
            plot.text(xmin-0.25, yvals[i], labels[i], fontsize=8, color='b', \
                horizontalalignment='right', verticalalignment='center')

    subplot.axis([xmin-1, xmax, -1, 2*N+1])
    subplot.axvline(0, linestyle='--', color='k', linewidth=1)

    subplot.axes.get_yaxis().set_visible(False)

    subplot.set_xlabel('proportion of QTLs', fontsize=10)
    for tick in subplot.get_xticklabels():
        tick.set_fontsize(10)

    return figure

def compute_proportions(filename):

    handle = open(filename, 'r')
    annot = cPickle.load(handle)
    test_statistics = cPickle.load(handle)
    variant_annotation = cPickle.load(handle)
    handle.close()

    prior_total = 0
    annot_prior = dict()
    pos_total = 0
    annot_posterior = dict()
    towrite = []
    for gene,value in test_statistics.iteritems():
        snps = value.keys()
        prior = np.array([value[snp][2] for snp in snps])
        prior_argmax = np.argmax(prior)
        prior_snp = snps[prior_argmax]
        for a in variant_annotation[prior_snp]:
            try:
                annot_prior[a] += value[prior_snp][2]*0.5
            except KeyError:
                annot_prior[a] = value[prior_snp][2]*0.5
            prior_total += value[prior_snp][2]*0.5
        posterior = np.array([value[snp][3] for snp in snps])
        pos_argmax = np.argmax(posterior)
        pos_snp = snps[pos_argmax]
        for a in variant_annotation[pos_snp]:
            try:
                annot_posterior[a] += value[pos_snp][3]*value[pos_snp][4]
            except KeyError:
                annot_posterior[a] = value[pos_snp][3]*value[pos_snp][4]
            pos_total += value[pos_snp][3]*value[pos_snp][4]
        towrite.append('%s\t%s\t%s\t%.6f'%(gene,','.join(snps),','.join(['%.6f'%value[snp][3] for snp in snps]),value[pos_snp][4]))

    labels = list(set(annot_prior.keys()).union(set(annot_posterior.keys())))
    pre_proportion = np.array([annot_prior[label] if annot_prior.has_key(label) else 0 for label in labels])/prior_total
    post_proportion = np.array([annot_posterior[label] if annot_posterior.has_key(label) else 0 for label in labels])/pos_total

    return pre_proportion, post_proportion, labels, towrite

def compute_weights(filename):

    handle = open(filename, 'r')
    annot = cPickle.load(handle)
    handle.close()
    weight = annot[1]
    stderr = annot[2]

    N = len(annot[0])
    annotdict = dict([(val,key) for key,val in annot[0].iteritems()])
    labels = [annotdict[i] for i in xrange(N)]

    return weight, stderr, labels

if __name__=="__main__":

    results_file = sys.argv[1]
    basename = os.path.splitext(results_file)[0]

    # get annotation weights    
    weight, stderr, labels = compute_weights(results_file)

    # plot annotation weights
    figure = plot_annotation_weights(weight, stderr, labels)
    annot_fig_file = '%s_annotation_weights.pdf'%basename
    figure.savefig(annot_fig_file, dpi=450)

    # write annotation weights
    weights_file = '%s_annotation_weights.txt'%basename
    handle = open(weights_file, 'w')
    handle.write("annotation\tweight\tstderr\n")
    for w,s,l in zip(weight,stderr,labels):
        handle.write("%s\t%.6f\t%.6f\n"%(l,w,s))
    handle.close()

    # compute proportions of QTLs in each annotation
    preprop, postprop, labels, towrite = compute_proportions(results_file)

    # write proportions
    proportion_file = '%s_probabilities.txt.gz'%basename
    handle = gzip.open(proportion_file, 'w')
    handle.write("locus\tsnps(comma-sep)\tsnp_posteriors(comma-sep)\tlocus_posterior\n")
    handle.write('\n'.join(towrite))
    handle.close()

    # plot proportions
    figure = plot_annotation_proportion(preprop, postprop, labels)
    prop_fig_file = '%s_annotation_proportions.pdf'%basename
    figure.savefig(prop_fig_file, dpi=450)
