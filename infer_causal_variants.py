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
    annot = [annotation.annot_labels, annotation.weights, annotation.stderr]

    # get per-gene and per-variant posteriors
    for datum,posterior in zip(data,posteriors):
        logBFmax = np.max(datum.logBF)
        prior = np.exp(datum.logBF-logBFmax)/np.sum(np.exp(datum.logBF-logBFmax))
        for snp,pri,pos in zip(datum.snps,prior,posterior.snp):
            test_statistics[datum.name][snp] = np.hstack((test_statistics[datum.name][snp],[pri,pos,posterior.gene]))

    # save data and results
    handle = open(options.output_file, 'w')
    cPickle.Pickler(handle,protocol=2).dump(annot)
    cPickle.Pickler(handle,protocol=2).dump(test_statistics)
    cPickle.Pickler(handle,protocol=2).dump(variant_annotations)
    handle.close()
