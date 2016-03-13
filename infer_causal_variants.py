import argparse
import cPickle
import warnings
import pdb
import sys

import numpy as np

import preprocess 
import bhmodel

# ignore warnings with these expressions
warnings.filterwarnings('ignore', '.*overflow encountered.*',)
warnings.filterwarnings('ignore', '.*divide by zero.*',)
warnings.filterwarnings('ignore', '.*invalid value.*',)

def parse_args():
    parser = argparse.ArgumentParser(description=" learns the parameters of riboHMM to infer translation "
                                     " from ribosome profiling data and RNA sequence data; "
                                    " RNA-seq data can also be used if available ")

    parser.add_argument("--mintol",
                        type=float,
                        default=1e-6,
                        help="convergence criterion for change in per-locus marginal likelihood (default: 1e-6)")

    parser.add_argument("--prior_var",
                        type=float,
                        default=1,
                        help="prior variance of effect sizes of causal variants (default: 1)")

    parser.add_argument("--output_file",
                        type=str,
                        default=None,
                        help="file to store the model parameters and per-variant and per-locus posteriors")

    parser.add_argument("--log_file",
                        type=str,
                        default=None,
                        help="file to store some statistics of the optimization algorithm")

    parser.add_argument("test_stat_file",
                        action="store",
                        help="file containing the test statistics")

    parser.add_argument("annot_files",
                        action="store",
                        nargs="+",
                        help="names of files containing the genomic and variant annotations")

    options = parser.parse_args()

    if options.output_file is None:
        raise NotImplementedError

    return options

def compute_posterior_enrichment(annotation, test_statistics, variant_annotation):

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

    labels = list(set(annot_prior.keys()).union(set(annot_posterior.keys())))
    pre_proportion = np.array([annot_prior[label] if annot_prior.has_key(label) else 0 for label in labels])/prior_total
    post_proportion = np.array([annot_posterior[label] if annot_posterior.has_key(label) else 0 for label in labels])/pos_total

    return pre_proportion, post_proportion, labels

if __name__=="__main__":

    options = parse_args()

    # load data
    test_statistics = preprocess.load_test_statistics(options.test_stat_file)
    print "loaded test statistics ..."
    genomic_annotations = preprocess.load_genomic_annotations(options.annot_files)
    print "loaded genomic annotations..."
    variant_annotations = preprocess.match_annotations_to_variants(test_statistics, genomic_annotations)
    print "preprocessed data ..."

    # learn model
    data, posteriors, annotation = bhmodel.learn_and_infer(test_statistics, variant_annotations, options.prior_var, options.mintol)

    # compute posterior enrichment
    preproportion, postproportion = compute_proportions(genomic_annotations, test_statistics, variant_annotations)

    # output annotation weights and standard errors,
    # along with their posterior enrichments
    annot_output_file = '%s_annotations.txt'%options.output_prefix
    handle = open(annot_output_file, 'w')
    handle.write("Annotation\tWeight\tStderr\tPosteriorEnrichment(woAnnotation)\tPosteriorEnrichment(withAnnotation)\n")
    for label,weight,stderr,pre,post in zip(annotation.annot_labels,annotation.weights,annotation.stderr,preproportion,postproportion):
        handle.write("%s\t%.8f\t%.8f\t%.8f\t%.8f\n"%(label,weight,stderr,pre,post))
    handle.close()

    # output per-locus QTL posterior
    locus_posterior_file = '%s_perlocus_posterior.txt.gz'%options.output_prefix
    handle = gzip.open(locus_posterior_file,'w')
    handle.write("Locus\tPosterior\n")
    for datum,posterior in zip(data,posteriors):
        handle.write("%s\t%.8f\n"%(datum.name,posterior.gene))
    handle.close()

    # output per-variant posterior of being the causal variant
    # with and without annotation information
    variant_posterior_file = '%s_pervariant_posterior.txt.gz'%options.output_prefix
    handle = gzip.open(variant_posterior_file, 'w')
    handle.write("Locus\tVariant\tPosterior(woAnnotation)\tPosterior(withAnnotation)\n")
    for datum,posterior in zip(data,posteriors):
        logBFmax = np.max(datum.logBF)
        prior = np.exp(datum.logBF-logBFmax)/np.sum(np.exp(datum.logBF-logBFmax))
        for snp,pri,pos in zip(datum.snps,prior,posterior.snp):
            handle.write("%s\t%s\t%.8f\t%.8f\n"%(datum.name, snp, pri, pos))
    handle.close()
