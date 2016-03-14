import argparse
import cPickle
import warnings
import gzip
import pdb
import os

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

    parser.add_argument("--output_prefix",
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

    if options.output_prefix is None:
        path = os.path.dirname(options.test_stat_file)
        stat_tag = os.path.basename(options.test_stat_file).split('.')[0]
        annot_tag = '_'.join([os.path.basename(filename).split('.')[0] for filename in options.annot_files])
        options.output_prefix = '%s/%s_%s_%s'%(path, stat_tag, annot_tag, str(options.prior_var))

    return options

def compute_posterior_enrichment(test_statistics, variant_annotation):

    annot_prior = dict()
    annot_posterior = dict()
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
        posterior = np.array([value[snp][3] for snp in snps])
        pos_argmax = np.argmax(posterior)
        pos_snp = snps[pos_argmax]
        for a in variant_annotation[pos_snp]:
            try:
                annot_posterior[a] += value[pos_snp][3]*value[pos_snp][4]
            except KeyError:
                annot_posterior[a] = value[pos_snp][3]*value[pos_snp][4]

    priortotal = np.sum(annot_prior.values())
    annot_prior = dict([(key,val/priortotal) for key,val in annot_prior.iteritems()])
    postotal = np.sum(annot_posterior.values())
    annot_posterior = dict([(key,val/postotal) for key,val in annot_posterior.iteritems()])

    return annot_prior, annot_posterior

if __name__=="__main__":

    options = parse_args()

    # load data
    test_statistics = preprocess.load_test_statistics(options.test_stat_file)
    print "%s\tloaded test statistics from %s"%(bhmodel.time(), options.test_stat_file)
    genomic_annotations = preprocess.load_genomic_annotations(options.annot_files)
    print "%s\tloaded genomic annotations from %s"%(bhmodel.time(), ",".join(options.annot_files))
    variant_annotations = preprocess.match_annotations_to_variants(test_statistics, genomic_annotations)
    print "%s\tpreprocessed data."%bhmodel.time()

    # learn model
    data, posteriors, annotation = bhmodel.learn_and_infer(test_statistics, variant_annotations, options.prior_var, options.mintol)

    # output locus QTL posterior
    locus_posterior_file = '%s_locus_posterior.txt.gz'%options.output_prefix
    handle = gzip.open(locus_posterior_file,'w')
    handle.write("Locus\tPosterior\n")
    for datum,posterior in zip(data,posteriors):
        handle.write("%s\t%.8f\n"%(datum.name,posterior.gene))
    handle.close()
    print "%s\toutput locus-level posteriors."%bhmodel.time()

    # output variant posterior of being the causal variant
    # with and without annotation information
    variant_posterior_file = '%s_variant_posterior.txt.gz'%options.output_prefix
    handle = gzip.open(variant_posterior_file, 'w')
    handle.write("Locus\tVariant\tPosterior(woAnnotation)\tPosterior(withAnnotation)\n")
    for datum,posterior in zip(data,posteriors):
        logBFmax = np.max(datum.logBF)
        prior = np.exp(datum.logBF-logBFmax)/np.sum(np.exp(datum.logBF-logBFmax))
        for snp,pri,pos in zip(datum.snps,prior,posterior.snp):
            handle.write("%s\t%s\t%.8f\t%.8f\n"%(datum.name, snp, pri, pos))
            test_statistics[datum.name][snp] = np.hstack((test_statistics[datum.name][snp], [pri, pos, posterior.gene]))
    handle.close()
    print "%s\toutput variant-level posteriors."%bhmodel.time()

    # compute posterior enrichment
    annot_prior, annot_posterior = compute_posterior_enrichment(test_statistics, variant_annotations)
    print "%s\tcomputed posterior enrichment within annotations."%bhmodel.time()

    # output annotation weights and standard errors,
    # along with their posterior enrichments
    annot_output_file = '%s_annotations.txt'%options.output_prefix
    handle = open(annot_output_file, 'w')
    handle.write("Annotation\tWeight\tStderr\tPosteriorEnrichment(woAnnotation)\tPosteriorEnrichment(withAnnotation)\n")
    for label,weight,stderr in zip(annotation.annot_labels,annotation.weights,annotation.stderr):
        try:
            prior = annot_prior[label]
        except KeyError:
            prior = 0
        try:
            posterior = annot_posterior[label]
        except KeyError:
            posterior = 0
        handle.write("%s\t%.8f\t%.8f\t%.8f\t%.8f\n"%(label,weight,stderr,prior,posterior))
    handle.close()
    print "%s\toutput posterior enrichment within annotations."%bhmodel.time()
    print "%s\tfinished successfully."%bhmodel.time()
