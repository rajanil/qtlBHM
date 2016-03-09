
import numpy as np
cimport numpy as np
from cpython cimport bool

cdef class Datum:

    cdef public long V
    cdef public str name
    cdef public list snps
    cdef public np.ndarray logBF

    cdef compute_bayes_factors(self, np.ndarray[np.float64_t, ndim=1] beta, np.ndarray[np.float64_t, ndim=1] stderr, prior_var)

cdef class Posterior:

    cdef public double gene
    cdef public np.ndarray snp

    cdef update(self, Datum datum, Annotation annotation)

cdef class Annotation:

    cdef public long N
    cdef public double log_prior_odds
    cdef public np.ndarray weights
    cdef public np.ndarray stderr
    cdef public dict annot_labels
    cdef public dict annotvalues

cdef tuple compute_func_grad(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors)

cdef tuple compute_func_grad_hess(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors)

cdef public double likelihood(list data, list posterior, Annotation annotation)
