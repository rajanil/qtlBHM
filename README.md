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
error messages, the compilation step should produce three files: 
bhmodel.c, bhmodel.h and bhmodel.so.

## Executing the code

The script you will need to execute is `infer_causal_variants.py`. To see command-line
options that need to be passed to the script, you can do the following:

    $ python infer_causal_variants.py

    runs qtlBHM, to infer causal variants driving molecular quantitative trait loci
    and identifies functional elements underlying these causal variants.
    qtlBHM requires summary statistics of association tests between the molecular phenotype
    for each locus and the genotypes of cis-variants, along with functional annotations 
    underlying each variant being tested.

    positional arguments:
      test_stat_file        gzipped text file containing summary statistics of
                            association tests between the molecular phenotypes of
                            a number of genetic loci and the genotype of variants
                            segregating in the sample. columns of the file should be
                            as follows: Locus_ID Variant Effect_size Standard_error
                            Variants should be specified as Chromosome.Position
                            (E.g., chr12.472368)
      annot_files           whitespace separated list of gzipped bed files 
                            containing genomic or variant annotations

    optional arguments:
      -h, --help            show this help message and exit
      --output_prefix OUTPUT_PREFIX
                            prefix for the output files
      --prior_var PRIOR_VAR prior variance of effect sizes of the causal
                            variants (default: 1.)
      --mintol MINTOL       convergence criterion for change in per-locus marginal
                            likelihood (default: 1e-6)
      --log_file LOG_FILE   file to store some statistics of the optimization algorithm

We will now describe in detail how to use this software using an example dataset of expression QTLs identified in human lymphoblastoid cell lines (LCLs), and chromHMM annotations for the tested variants. The test statistics and annotations are provided in `test/`.

### Key Inputs

The key inputs that need to be passed to this script are
+   the gzipped file containing the test statistics (effect sizes and standard errors)
+   the gzipped file(s) containing genomic and/or variant annotations

    *Note: these inputs are positional arguments and the files must be specified in the correct order (as shown above).*

The gzipped file of test statistics should have the following format 
(see `test/statistics.txt.gz` for an example; the header should not be included in the file).

    Locus               Variant         EffectSize  StandardError
    ENSG00000183020.9   chr22.19910762  0.00182715  0.142038689611108
    ENSG00000110344.5   chr22.45258495  -0.675393   0.226839423859863

The gzipped annotation file should have the following format (the header should not be included in the file).
(see `test/annotation.bed.gz` for an example; the header should not be included in the file).

    Chromosome  Start   End         AnnotationLabel
    chr22   51221735    51223534    Active_promoter_LCL
    chr22   17639200    17639399    Weak_promoter_LCL

In the above formats, positions can be 0-based or 1-based, as long as they are the same in both files.

### Learning and Inference

Learning and inference can be performed by passing the following arguments (using test data as examples).

    python infer_causal_variants.py --output_prefix test/results test/statistics.txt.gz test/annotations.bed.gz

This run will output three files: 
  1. `test/results_locus_posterior.txt.gz` -- this file lists the probability that each tested locus is a QTL
  2. `test/results_variant_posterior.txt.gz` -- this file lists the probability that each variant is the causal variant for a specific locus.A probability is computed with and without prior weighting of variants using annotations.
  3. `test/results_annotations.txt` -- this file lists the weights and standard errors for each annotation, and the posterior enrichment for each annotation. A posterior enrichment is computed with and without prior weighting of variants using annotations.
