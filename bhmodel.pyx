import numpy as np
cimport numpy as np
from math import log, exp
import random
import cvxopt as cvx
from cvxopt import solvers
import datetime
import pdb

solvers.options['maxiters'] = 20
solvers.options['show_progress'] = False

time = lambda: datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

cdef double EPS = np.finfo(np.double).tiny

cdef np.ndarray nplogvec(np.ndarray x):

    x[x<EPS] = EPS
    return np.log(x)

cdef double nplog(double x):

    return np.log(max([x,EPS]))

def write_log(log_handle, towrite):

    log_handle.write(towrite+'\n')
    print towrite

cdef class Datum:

    def __cinit__(self, str gene, dict associations, Prior prior):

        cdef str key
        self.name = gene
        self.snps = associations.keys()
        self.V = len(self.snps)
        self.beta = np.array([associations[key][0] for key in self.snps])
        self.stderr = np.array([associations[key][1] for key in self.snps])
        self.logBF = np.empty((self.V,), dtype=np.float64)
        self.compute_bayes_factors(prior)

    cdef compute_bayes_factors(self, Prior prior):
        """
        compute an approximate log Bayes Factor, as detailed
        in Wakefield 2008
        We use -log(ABF) described in fomula (2) in Wakefield 2008.
        We flip the sign so that positive values indicate support for the
        alternative rather than null.
        """

        cdef np.ndarray zscore, r

        zscore = self.beta/self.stderr
        r = prior.var / (self.stderr**2 + prior.var)
        self.logBF = 0.5*nplogvec(1-r) + 0.5*r*zscore**2

cdef class Posterior:

    def __cinit__(self, Datum datum):

        self.gene = 0
        self.snp = np.empty((datum.V,),dtype=float)

    cdef update(self, Datum datum, Annotation annotation, Prior prior):

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
            chi[v] = datum.logBF[v] + chi[v] - nplog(chisum)
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
        self.gene = prior.gene_prior_logodds + np.sum(self.snp*(chi-nplogvec(self.snp)))
        try:
            self.gene = 1./(1+exp(-1*self.gene))
        except OverflowError:
            self.gene = 0
        if self.gene<0 or self.gene>1:
            print self.gene
            pdb.set_trace()

cdef class Prior:

    def __cinit__(self):

        # prior probability that a gene has a causal QTL = 0.1
        self.gene_prior_logodds = -1*nplog(9)

        # initial variance of causal QTLs is drawn from standard
        # inverse-gamma distribution
        self.var = 1./np.random.gamma(1,1)

    def update(self, list data, list posteriors):

        #self.update_gene_prior(posteriors)

        self.update_var(data, posteriors)

    def update_var(self, list data, list posteriors):

        x_init = np.array([self.var])
        self.var = optimize_prior_variance(x_init, data, posteriors)

    cdef update_gene_prior(self, list posteriors):

        # update prior probability that a gene has a causal eQTL
        self.gene_prior_logodds = np.sum([posterior.gene for posterior in posteriors])/len(posteriors)
        self.gene_prior_logodds = nplog(self.gene_prior_logodds/(1-self.gene_prior_logodds))

def optimize_prior_variance(x_init, data, posteriors):

    # function that computes likelihood, gradient and hessian
    def F(x=None, z=None):

        if x is None:
            return 0, cvx.matrix(x_init)

        xx = np.array(x).ravel().astype(np.float64)

        if z is None:
            results = compute_func_grad_priorvar(xx, data, posteriors)
            f = results[0]
            Df = results[1]
            return cvx.matrix(f), cvx.matrix(Df)
        else:
            results = compute_func_grad_hess_priorvar(xx, data, posteriors)
            f = results[0]
            Df = results[1]
            Hf = z[0]*results[2]
            return cvx.matrix(f), cvx.matrix(Df), cvx.matrix(Hf)

    G = cvx.matrix(np.array([[-1.]]))
    h = cvx.matrix(np.array([[0.]]))
    try:
        solution = solvers.cp(F, G=G, h=h)
        x_final = np.array(solution['x'])[0]
    except ValueError:
        x_final = x_init[0]
    return x_final

cdef tuple compute_func_grad_priorvar(np.ndarray[np.float64_t, ndim=1] xx, list data, list posteriors):

    cdef double f
    cdef np.ndarray[np.float64_t, ndim=2] df
    cdef Posterior posterior
    cdef Datum datum

    f = 0.
    df = np.zeros((1,1), dtype='float')

    for datum,posterior in zip(data,posteriors):

        f = f + posterior.gene * np.sum(posterior.snp * \
            (2*nplogvec(datum.stderr) - nplogvec(xx+datum.stderr**2) + xx/(xx+datum.stderr**2) * (datum.beta/datum.stderr)**2))

        df[0,0] = df[0,0] + posterior.gene * np.sum(posterior.snp * \
                ((datum.beta**2 - datum.stderr**2 - xx)/(xx + datum.stderr**2)**2))

    return -1*f, -1*df

cdef tuple compute_func_grad_hess_priorvar(np.ndarray[np.float64_t, ndim=1] xx, list data, list posteriors):

    cdef double f
    cdef np.ndarray[np.float64_t, ndim=2] df, hess
    cdef Posterior posterior
    cdef Datum datum

    f = 0.
    df = np.zeros((1,1), dtype='float')
    hess = np.zeros((1,1), dtype='float')

    for datum,posterior in zip(data,posteriors):

        f = f + posterior.gene * np.sum(posterior.snp * \
            (2*nplogvec(datum.stderr) - nplogvec(xx+datum.stderr**2) + xx/(xx+datum.stderr**2) * (datum.beta/datum.stderr)**2))

        df[0,0] = df[0,0] + posterior.gene * np.sum(posterior.snp * \
                  ((datum.beta**2 - datum.stderr**2 - xx)/(xx + datum.stderr**2)**2))

        hess[0,0] = hess[0,0] + posterior.gene * np.sum(posterior.snp * \
                    (xx + datum.stderr**2 - 2*datum.beta**2)/(xx + datum.stderr**2)**3)

    return -1*f, -1*df, -1*hess

cdef class Annotation:

    def __cinit__(self, dict annot_values):

        cdef long i
        cdef str snp, label, v
        cdef list val, annot_labels

        # only select annotations that cover a minimum number of variants (100)
        all_labels = [v for val in annot_values.values() for v in val]
        annot_labels = list(set(all_labels))
        annot_labels.sort()

        self.N = len(annot_labels)
        self.weights = np.zeros((self.N,),dtype=float)
        self.stderr = np.zeros((self.N,),dtype=float)
        self.annot_labels = dict([(label,i) for i,label in enumerate(annot_labels)])

        self.annotvalues = dict()
        for snp,val in annot_values.iteritems():
            self.annotvalues[snp] = [self.annot_labels[v] for v in val 
                                     if self.annot_labels.has_key(v)]
            self.annotvalues[snp].sort()

    def update(self, list data, list posteriors):

        # run solver to update annotation weights
        x_init = self.weights.reshape(self.N,1)
        self.weights = optimize_annotation_weights(x_init, data, self.annotvalues, posteriors)

    def compute_stderr(self, list data, list posteriors):

        results = compute_func_grad_hess_weights(self.weights, data, self.annotvalues, posteriors)
        self.stderr = np.sqrt(np.diag(np.linalg.inv(results[2])))


def optimize_annotation_weights(x_init, data, annotvalues, posteriors):

    # function that computes likelihood, gradient and hessian
    def F(x=None, z=None):

        if x is None:
            return 0, cvx.matrix(x_init)

        xx = np.array(x).ravel().astype(np.float64)

        if z is None:
            results = compute_func_grad_weights(xx, data, annotvalues, posteriors)
            f = results[0]
            Df = results[1]
            return cvx.matrix(f), cvx.matrix(Df)
        else:
            results = compute_func_grad_hess_weights(xx, data, annotvalues, posteriors)
            f = results[0]
            Df = results[1]
            Hf = z[0]*results[2]
            return cvx.matrix(f), cvx.matrix(Df), cvx.matrix(Hf)

    solution = solvers.cp(F)
    x_final = np.array(solution['x']).ravel()
    return x_final

cdef tuple compute_func_grad_weights(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors):

    cdef long a, v, V, index
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

        f = f + posterior.gene * (total - nplog(chisum) - chimax)
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

cdef tuple compute_func_grad_hess_weights(np.ndarray[np.float64_t, ndim=1] xx, list data, dict annotvalues, list posteriors):

    cdef long a, b, v, V, aidx, bidx, I, J, i, j, index
    cdef double f, chimax, chisum, total, u, w
    cdef list keys, annotval
    cdef dict expect
    cdef np.ndarray[np.float64_t, ndim=1] df, chi
    cdef np.ndarray[np.float64_t, ndim=2] hess
    cdef Posterior posterior
    cdef Datum datum

    V = xx.size
    J = len(posteriors)
    f = 0.
    df = np.zeros((V,), dtype='float')
    hess = 1./J*np.eye(V)

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

        f = f + posterior.gene * (total - nplog(chisum) - chimax)
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
                
cdef double likelihood(list data, list posteriors, Annotation annotation, Prior prior):

    cdef v
    cdef int N
    cdef double L, chisum, temp
    cdef np.ndarray chi
    cdef Datum datum
    cdef Posterior posterior

    L = 0.
    N = len(data)
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

        chisum = nplog(chisum)
        temp = 0
        for v from 0 <= v < datum.V:
            chi[v] = datum.logBF[v] + nplog(chi[v]) - chisum
            temp = temp + posterior.snp[v] * chi[v]
            try:
                temp = temp - posterior.snp[v]*nplog(posterior.snp[v])
            except ValueError:
                pass

        L = L + posterior.gene * (temp + prior.gene_prior_logodds)
        L = L - posterior.gene*nplog(posterior.gene) - (1-posterior.gene)*nplog(1-posterior.gene)

    return L/N

def learn_and_infer(dataQTL, snp_annotation, log_handle, reltol):

    cdef list data, posteriors
    cdef Annotation annotation
    cdef Datum datum
    cdef Prior prior
    cdef Posterior posterior

    # initialize variables and prior
    prior = Prior()
    annotation = Annotation(snp_annotation)

    # initialize data
    data = [Datum(key,value,prior) for key,value in dataQTL.iteritems()]

    # initialize latent variables and likelihood
    posteriors = [Posterior(datum) for datum in data]
    [posterior.update(datum, annotation, prior) for datum,posterior in zip(data,posteriors)]
    L = likelihood(data, posteriors, annotation, prior)
    write_log(log_handle, "%s\t%d genes; %d variants; %d annotations"%(time(),len(data),len(annotation.annotvalues),annotation.N))
    write_log(log_handle, "%s\tInitial likelihood = %.6e"%(time(), L))
    dL = np.inf

    # infer causal variants
    while np.abs(dL)>reltol:

        # update annotation parameters 
        annotation.update(data, posteriors)

        # update prior
        prior.update(data, posteriors)

        # update Bayes Factors and posterior
        [datum.compute_bayes_factors(prior) for datum in data]
        [posterior.update(datum, annotation, prior) for datum,posterior in zip(data,posteriors)]

        # compute likelihood
        newL = likelihood(data, posteriors, annotation, prior)

        dL = (newL-L)/np.abs(L)
        L = newL
        write_log(log_handle, "%s\tLikelihood = %.6e; Relative change = %.6e"%(time(),L,dL))

    annotation.compute_stderr(data, posteriors)
    return data, posteriors, annotation, prior
