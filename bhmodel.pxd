
import numpy as np
cimport numpy as np
from cpython cimport bool

cdef class Datum:

    cdef public long V
    cdef public str name
    cdef public list snps
    cdef public np.ndarray beta, stderr, logBF

    cdef compute_bayes_factors(self, Prior prior)

cdef class Posterior:

    cdef public double gene
    cdef public np.ndarray snp

    cdef update(self, Datum datum, Annotation annotation, Prior prior)

cdef class Prior:

    cdef public double var
    cdef public double gene_prior_logodds

    cdef update_gene_prior(self, list posteriors)

cdef class Annotation:

    cdef public long N
    cdef public double log_prior_odds
    cdef public np.ndarray weights
    cdef public np.ndarray stderr
    cdef public dict annot_labels
    cdef public dict annotvalues

cdef tuple compute_func_grad_weights(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors)

cdef tuple compute_func_grad_hess_weights(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors)

cdef tuple compute_func_grad_priorvar(np.ndarray[np.float64_t, ndim=1] xx, list data, list posteriors)

cdef tuple compute_func_grad_hess_priorvar(np.ndarray[np.float64_t, ndim=1] xx, list data, list posteriors)

cdef public double likelihood(list data, list posterior, Annotation annotation, Prior prior)

cdef public double nplog(double x)

cdef public np.ndarray nplogvec(np.ndarray x)
