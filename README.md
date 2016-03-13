# qtlBHM

**qtlBHM** is an algorithm for identifying functional elements underlying putatively causal variants
associated with any molecular trait (gene expression, chromatin accessibility, histone modifications, etc) 
and, simultaneously, resolving the causal variant from a set of associated variants. The algorithm is
written in Python2.x.

The [Bayesian hierarchical model underlying qtlBHM]() identifies causal variants driving a
quantitative trait locus (QTL), from a set of variants associated with the locus, 
by sharing information across different loci and quantifying the
enrichment of strongly associated variants within different functional elements. The model starts from
summary statistics (effect size and standard error) of association across individuals 
between the molecular phenotype of a locus and a set of segregating variants, and
the functional annotations underlying these variants. Inference under the model gives
the posterior probability that the locus is a QTL and the posterior probability that 
each of the tested variants is the causal variant driving that QTL.

This repo contains set of scripts to load the data, run the algorithm, and visualize the results. 
This document summarizes how to download and setup this software package and provides instructions 
on how to run the software on a test dataset of variants associated with gene expression and
their underlying functional annotations.

## Dependencies

qtlBHM depends on
+ [Numpy](http://www.numpy.org/)
+ [Cvxopt](http://www.cvxopt.org/)
+ [Matplotlib](http://matplotlib.org/)

A number of python distributions already have Numpy and Matplotlib packaged in them. It is also
straightforward to install all these dependencies
 (1) using package managers for MACOSX and several Linux distributions,
 (2) from platform-specific binary packages, and
 (3) directly from source

## Getting the source code

To obtain the source code from github, let us assume you want to clone this repo into a
directory named `proj`:

    mkdir ~/proj
    cd ~/proj
    git clone https://github.com/rajanil/qtlBHM

To retrieve the latest code updates, you can do the following:

    cd ~/proj/qtlBHM
    git fetch
    git merge origin/master

Since the algorithm is written in Cython, qtlBHM will have to be compiled into a shared object in 
the following way, before it can be executed.

    python setup.py build_ext --inplace

This step will generate a number of warning messages that can be ignored. If there are no 
error messages, the compilation step should produce two files bhmodel.c and bhmodel.so.


