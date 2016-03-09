import gzip
import pdb

import numpy as np

def load_test_statistics(filename):

    teststats = dict()

    # open file 
    if filename.lower().endswith(('.gz','.zip')):
        handle = gzip.open(filename, 'r')
    else:
        handle = open(filename, 'r')

    # read in effect sizes and standard errors
    for line in handle:
        row = line.strip().split()
        try:
            teststats[row[0]][row[1]] = np.array([float(row[2]), float(row[3])])
        except KeyError:
            teststats[row[0]] = {row[1]: np.array([float(row[2]), float(row[3])])}

    handle.close()

    return teststats

def load_genomic_annotations(filename):

    annotations = dict()

    # open file
    if filename.lower().endswith(('.gz','.zip')):
        handle = gzip.open(filename, 'r')
    else:
        handle = open(filename, 'r')

    # read in annotation coordinates and labels
    for line in handle:
        row = line.strip().split()
        # check if annotation exists in dictionary
        try:
            annotations[row[3]]
        except KeyError:
            annotations[row[3]] = dict()
        # add entry
        try:
            annotations[row[3]][row[0]].append([int(row[1]), int(row[2])])
        except KeyError:
            annotations[row[3]][row[0]] = [[int(row[1]), int(row[2])]]

    handle.close()

    # order annotations by location
    for label in annotations.keys():
        for chrom in annotations[label].keys():
            values = np.array(annotations[label][chrom])
            order = np.argsort(values[:,0])
            annotations[label][chrom] = values[order,:]

    return annotations

def match_annotations_to_variants(test_statistics, genomic_annotations):

    variant_annotations = dict()
    for locus,stats in test_statistics.iteritems():
        for variant,varstat in stats.iteritems():
            chrom, position = variant.split('.')
            position = int(position)

            try:
                variant_annotations[variant]
            except KeyError:
                variant_annotations[variant] = [label for label,values in genomic_annotations.iteritems() 
                                                if np.any(np.logical_and(values[chrom][:,0]<=position,values[chrom][:,1]>position))]

    return variant_annotations                
