import sys
import pdb

import numpy as np
import matplotlib.pyplot as plot

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

    order = np.argsort(postprop)
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

def get_annotation_weights_proportions(filename):

    handle = open(filename, 'r')
    header = handle.next().strip().split()
    data = [line.strip().split() for line in handle]
    handle.close()
    labels = [dat[0] for dat in data]
    weight = np.array([dat[1] for dat in data]).astype('float')
    stderr = np.array([dat[2] for dat in data]).astype('float')
    preprop = np.array([dat[3] for dat in data]).astype('float')
    postprop = np.array([dat[4] for dat in data]).astype('float')

    return weight, stderr, preprop, postprop, labels

if __name__=="__main__":

    results_prefix = sys.argv[1]

    # get annotation weights    
    weight, stderr, preprop, postprop, labels = get_annotation_weights_proportions(results_prefix+"_annotations.txt")

    # plot annotation weights
    figure = plot_annotation_weights(weight, stderr, labels)
    annot_fig_file = '%s_annotation_weights.pdf'%results_prefix
    figure.savefig(annot_fig_file, dpi=450)

    # plot proportions
    figure = plot_annotation_proportion(preprop, postprop, labels)
    prop_fig_file = '%s_annotation_proportions.pdf'%results_prefix
    figure.savefig(prop_fig_file, dpi=450)
