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
                        help="prior variance on effect sizes (default: 1)")

    parser.add_argument("--output_file",
                        type=str,
                        default=None,
                        help="file to store the model parameters and per-variant posteriors")

    parser.add_argument("test_stat_file",
                        action="store",
                        help="file containing the test statistics")

    parser.add_argument("annot_file",
                        action="store",
                        help="file containing the genomic annotations")

    options = parser.parse_args()

    if options.output_file is None:
        raise NotImplementedError

    return options

if __name__=="__main__":

    options = parse_args()

    # load data
    test_statistics = preprocess.load_test_statistics(options.test_stat_file)
    genomic_annotations = preprocess.load_genomic_annotations(options.annot_file)
    variant_annotations = preprocess.match_annotations_to_variants(test_statistics, genomic_annotations)
    pdb.set_trace()

    # learn model
    data, posteriors, annotation = bhmodel.learn_and_infer(test_statistics, variant_annotations, options.prior_var, options.reltol)
    annot = [annotation.annot_labels, annotation.weights, annotation.stderr]

    # get per-gene and per-variant posteriors
    for datum,posterior in zip(data,posteriors):
        for snp,pos in zip(datum.snps,posterior.snp):
            rawdata[datum.name][snp] = np.hstack((rawdata[datum.name][snp],[pos,posterior.gene]))

    # save data and results
    handle = open(option.output_file, 'w')
    cPickle.Pickler(handle,protocol=2).dump(annot)
    cPickle.Pickler(handle,protocol=2).dump(rawdata)
    handle.close()
