import gzip
import pdb
import time

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

def load_genomic_annotations(filenames):

    annotations = dict()

    for filename in filenames:

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

    # get all variants
    variants = dict()
    for locus,stats in test_statistics.iteritems():
        for variant,varstat in stats.iteritems():
            chrom, position = variant.split('.')
            position = int(position)
            try:
                variants[chrom].append(position)
            except KeyError:
                variants[chrom] = [position]

    # initialize annotation dictionary
    variant_annotations = dict()
    for chrom in variants.keys():
        variants[chrom] = np.unique(variants[chrom])
        variant_annotations.update(dict([('%s.%d'%(chrom,pos),[]) for pos in variants[chrom]]))

    for label,values in genomic_annotations.iteritems():
        for chrom in values.keys():
            try:
                variants[chrom]
            except KeyError:
                continue
            chrmax = max([values[chrom].max(), variants[chrom].max()])+1
            mask = np.zeros((chrmax,), dtype='bool')
            [mask.__setslice__(val[0], val[1], True) for val in values[chrom]]
            vars = ['%s.%d'%(chrom,pos) for pos in variants[chrom][mask[variants[chrom]]]]
            ig = [variant_annotations[var].append(label) for var in vars]

    return variant_annotations
