import numpy as np
cimport numpy as np
from math import log, exp
import random
import cvxopt as cvx
from cvxopt import solvers
import time
import pdb

solvers.options['maxiters'] = 20
solvers.options['show_progress'] = False

cdef class Datum:

    def __cinit__(self, gene, associations, prior_var):

        cdef str key
        cdef np.ndarray[np.float64_t, ndim=1] beta, stderr
        self.name = gene
        self.snps = associations.keys()
        self.V = len(self.snps)
        beta = np.array([associations[key][1] for key in self.snps])
        stderr = np.array([associations[key][2] for key in self.snps])
        self.logBF = np.empty((self.V,), dtype=np.float64)
        self.compute_bayes_factors(beta, stderr, prior_var)

    cdef compute_bayes_factors(self, np.ndarray[np.float64_t, ndim=1] beta, \
    np.ndarray[np.float64_t, ndim=1] stderr, prior_var):
        """
        compute an approximate log Bayes Factor, as detailed
        in Wakefield 2008
        We use -log(ABF) described in fomula (2) in Wakefield 2008.
        We flip the sign so that positive values indicate support for the
        alternative rather than null.
        """

        cdef long v
        cdef double zscore, r

        for v from 0 <= v < self.V:

            zscore = beta[v]/stderr[v]
            r = prior_var / (stderr[v]**2 + prior_var)
            self.logBF[v] = 0.5*log(1-r) + 0.5*r*zscore**2

cdef class Posterior:

    def __cinit__(self, Datum datum):

        self.gene = 0.1*random.random()
        self.snp = np.empty((datum.V,),dtype=float)

    cdef update(self, Datum datum, Annotation annotation):

        cdef long v, idx
        cdef double snpmax, chimax, chisum, snpsum
        cdef np.ndarray chi

        # update SNP posteriors
        chi = np.zeros((datum.V,), dtype=np.float64)
        chimax = -np.inf
        for v from 0 <= v < datum.V:
            for idx in annotation.annotvalues[datum.snps[v]]:
                chi[v] = chi[v] + annotation.weights[idx]
            chimax = max([chimax, chi[v]])

        chisum = 0
        for v from 0 <= v < datum.V:
            chi[v] = chi[v]-chimax
            chisum = chisum + exp(chi[v])

        snpmax = -np.inf
        for v from 0 <= v < datum.V:
            chi[v] = datum.logBF[v] + chi[v] - log(chisum)
            snpmax = max([snpmax, chi[v]])

        snpsum = 0
        for v from 0 <= v < datum.V:
            self.snp[v] = exp(chi[v]-snpmax)
            snpsum = snpsum + self.snp[v]
        self.snp = self.snp/snpsum
        if np.any(self.snp<0) or np.any(self.snp>1):
            print self.snp
            pdb.set_trace()

        # update gene posteriors
        self.gene = annotation.log_prior_odds + np.sum(self.snp*(chi-np.nan_to_num(np.log(self.snp))))
        try:
            self.gene = 1./(1+exp(-1*self.gene))
        except OverflowError:
            self.gene = 0
        if self.gene<0 or self.gene>1:
            print self.gene
            pdb.set_trace()

cdef class Annotation:

    def __cinit__(self, dict annot_values):

        cdef long i
        cdef str snp, label, v
        cdef list val, annot_labels

        self.log_prior_odds = random.random()
        self.log_prior_odds = log(self.log_prior_odds/(1-self.log_prior_odds))

        annot_labels = list(set([v for val in annot_values.values() for v in val]))
        annot_labels.sort()
        self.N = len(annot_labels)
        self.weights = np.zeros((self.N,),dtype=float)
        self.stderr = np.zeros((self.N,),dtype=float)
        self.annot_labels = dict([(label,i) for i,label in enumerate(annot_labels)])

        self.annotvalues = dict()
        for snp,val in annot_values.iteritems():
            self.annotvalues[snp] = [self.annot_labels[v] for v in val]
            self.annotvalues[snp].sort()

    def update(self, list data, list posteriors):

        # select subset of genes to update weights
        #subset = [(d,p) for d,p in zip(data,posteriors) if random.random()<0.1]
        #data = [s[0] for s in subset]
        #posteriors = [s[1] for s in subset]

        # run solver
        x_init = self.weights.reshape(self.N,1)
        self.weights = optimize_annotation_weights(x_init, data, self.annotvalues, posteriors)

    def compute_stderr(self, list data, list posteriors):

        results = compute_func_grad_hess(self.weights, data, self.annotvalues, posteriors)
        self.stderr = np.sqrt(np.diag(np.linalg.inv(results[2])))


def optimize_annotation_weights(x_init, data, annotvalues, posteriors):

    # function that computes likelihood, gradient and hessian
    def F(x=None, z=None):

        if x is None:
            return 0, cvx.matrix(x_init)

        xx = np.array(x).ravel().astype(np.float64)

        if z is None:
            results = compute_func_grad(xx, data, annotvalues, posteriors)
            f = results[0]
            Df = results[1]
            return cvx.matrix(f), cvx.matrix(Df)
        else:
            results = compute_func_grad_hess(xx, data, annotvalues, posteriors)
            f = results[0]
            Df = results[1]
            Hf = z[0]*results[2]
            return cvx.matrix(f), cvx.matrix(Df), cvx.matrix(Hf)

    solution = solvers.cp(F)
    x_final = np.array(solution['x']).ravel()
    return x_final

cdef tuple compute_func_grad(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors):

    cdef long a, v, V
    cdef double f, chimax, chisum, total, u
    cdef dict expect
    cdef np.ndarray[np.float64_t, ndim=1] df, chi
    cdef Posterior posterior
    cdef Datum datum

    V = xx.size
    f = 0.
    df = np.zeros((V,), dtype='float')
    for datum,posterior in zip(data,posteriors):

        # update SNP posteriors
        total = 0
        expect = dict()
        chi = np.zeros((datum.V,), dtype=np.float64)
        chimax = -np.inf
        for v from 0 <= v < datum.V:
            for a in annotvalues[datum.snps[v]]:
                chi[v] = chi[v] + xx[a]
            total = total + chi[v]*posterior.snp[v]
            chimax = max([chimax, chi[v]])

        chisum = 0
        for v from 0 <= v < datum.V:
            chi[v] = exp(chi[v]-chimax)
            chisum = chisum + chi[v]
        chi = chi/chisum

        f = f + posterior.gene * (total - log(chisum) - chimax)
        for v from 0 <= v < datum.V:
            u = posterior.gene * posterior.snp[v]
            for a in annotvalues[datum.snps[v]]:
                df[a] = df[a] + u
                try:
                    expect[a] = expect[a] + chi[v]
                except KeyError:
                    expect[a] = chi[v]
        
        for a,u in expect.iteritems():
            df[a] = df[a] - posterior.gene * u

    return -1*f, -1*df.reshape(1,V)

cdef tuple compute_func_grad_hess(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors):

    cdef long a, b, v, V, aidx, bidx, I, J, i, j
    cdef double f, chimax, chisum, total, u, w
    cdef list keys, annotval
    cdef dict expect
    cdef np.ndarray[np.float64_t, ndim=1] df, chi
    cdef np.ndarray[np.float64_t, ndim=2] hess
    cdef Posterior posterior
    cdef Datum datum

    V = xx.size
    f = 0.
    df = np.zeros((V,), dtype='float')
    hess = 0.01*np.eye(V)
    for datum,posterior in zip(data,posteriors):

        # update SNP posteriors
        total = 0
        chi = np.zeros((datum.V,), dtype=np.float64)
        chimax = -np.inf
        for v from 0 <= v < datum.V:
            for a in annotvalues[datum.snps[v]]:
                chi[v] = chi[v] + xx[a]
            total = total + chi[v]*posterior.snp[v]
            chimax = max([chimax, chi[v]])

        chisum = 0
        for v from 0 <= v < datum.V:
            chi[v] = exp(chi[v]-chimax)
            chisum = chisum + chi[v]
        chi = chi/chisum

        f = f + posterior.gene * (total - log(chisum) - chimax)
        expect = dict()
        for v from 0 <= v < datum.V:
            u = posterior.gene * posterior.snp[v]
            w = posterior.gene * chi[v]
            annotval = annotvalues[datum.snps[v]]
            J = len(annotval)
            for i from 0 <= i < J:
                a = annotval[i]
                df[a] = df[a] + u
                try:
                    expect[a] = expect[a] + chi[v]
                except KeyError:
                    expect[a] = chi[v]
                hess[a,a] = hess[a,a] + w
                for j from i < j < J:
                    b = annotval[j]
                    hess[a,b] = hess[a,b] + w

        keys = expect.keys()
        J = len(keys)
        for i from 0 <= i < J:
            a = keys[i]
            df[a] = df[a] - posterior.gene * expect[a]
            hess[a,a] = hess[a,a] - posterior.gene * expect[a]**2
            for j from i < j < J:
                b = keys[j]
                hess[a,b] = hess[a,b] - posterior.gene * expect[a] * expect[b]

    hess = hess + hess.T
    hess[range(V),range(V)] = 0.5*np.diag(hess)

    return -1*f, -1*df.reshape(1,V), hess
                
cdef double likelihood(list data, list posteriors, Annotation annotation):

    cdef v
    cdef double L, prior, chisum, temp
    cdef np.ndarray chi
    cdef Datum datum
    cdef Posterior posterior

    L = 0.
    for datum,posterior in zip(data,posteriors):

        # update SNP posteriors
        chi = np.zeros((datum.V,), dtype=np.float64)
        chimax = -np.inf
        for v from 0 <= v < datum.V:
            for val in annotation.annotvalues[datum.snps[v]]:
                chi[v] = chi[v] + annotation.weights[val]
            chimax = max([chimax, chi[v]])

        chisum = 0
        for v from 0 <= v < datum.V:
            chi[v] = exp(chi[v]-chimax)
            chisum = chisum + chi[v]

        chisum = log(chisum)
        temp = 0
        for v from 0 <= v < datum.V:
            chi[v] = datum.logBF[v] + log(chi[v]) - chisum
            temp = temp + posterior.snp[v] * chi[v]
            try:
                temp = temp - posterior.snp[v]*log(posterior.snp[v])
            except ValueError:
                pass

        L = L + posterior.gene * (temp + annotation.log_prior_odds)
        try:
            L = L - posterior.gene*log(posterior.gene) - (1-posterior.gene)*log(1-posterior.gene)
        except ValueError:
            pass

    return L

def learn_and_infer(dataQTL, snp_annotation, prior_var, reltol):

    cdef list data, posteriors
    cdef Annotation annotation
    cdef Datum datum
    cdef Posterior posterior

    # initialize data and variables
    data = [Datum(key,value,prior_var) for key,value in dataQTL.iteritems()]
    annotation = Annotation(snp_annotation)
    posteriors = [Posterior(datum) for datum in data]
    for datum,posterior in zip(data,posteriors):
        posterior.update(datum, annotation)
    L = likelihood(data, posteriors, annotation)
    print "%d genes"%len(data)
    print "%d variants"%(len(annotation.annotvalues))
    print "%d annotations"%annotation.N
    print "Initial likelihood = %.3e"%L
    dL = np.inf

    # infer causal variants
    while np.abs(dL)>reltol:

        starttime = time.time()

        # update posterior
        starttime = time.time()
        for datum,posterior in zip(data,posteriors):
            posterior.update(datum, annotation)
        print "updated posteriors in %.2f secs"%(time.time()-starttime)

        # update parameters
        starttime = time.time()
        annotation.update(data, posteriors)
        print "updated parameters in %.2f secs"%(time.time()-starttime)

        # compute likelihood
        starttime = time.time()
        newL = likelihood(data, posteriors, annotation)
        print "computed likelihood in %.2f secs"%(time.time()-starttime)

        dL = (newL-L)/np.abs(L)
        L = newL
        print "Likelihood = %.3e"%L, "Relative change = %.3e"%dL

    annotation.compute_stderr(data, posteriors)
    return data, posteriors, annotation
