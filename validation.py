#!/usr/bin/env python

from functools import wraps
from time import perf_counter
from typing import Callable

from hypothesis import given, settings
import hypothesis.strategies as st

import random
from itertools import product, permutations, combinations

def timing(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        start = perf_counter()
        result = f(*args, **kwargs)
        end = perf_counter()
        print(f'{f.__name__: <20}: {end-start:.2E}')
        return result
    return wrapper


def gen_occ(n_alpha,n_beta, seed = None, sparsity=0.20):
    random.seed(seed)
   
    if n_alpha < n_beta:
        n_alpha, n_beta = n_beta, n_alpha

    n_ref = n_alpha + n_beta
    l = []
    while True:
       n = sum(l)
       if n  == n_ref:
          break
       elif n + 1 == n_ref:
          o = 1
       elif sum(i for i in l if i == 2) >= min(n_alpha, n_beta):
          o = 1
       else: 
          o = random.randint(1,2)

       l.append(o)

    l += [0]*random.randint(0,int(len(l)*sparsity))
    random.shuffle(l) 
    return n_alpha, n_beta, l

@timing
def naive(occ, n_alpha, n_beta):

    n_orbital = len(occ)
    p =  list(product([0,1],repeat=n_orbital))
    # Filter alpha
    p_alpha = filter(lambda o: sum(o) == n_alpha, p )
    p_beta =  filter(lambda o: sum(o) == n_beta, p )

    # All combinaison
    q = product(p_alpha, p_beta)

    # No filter more than double occupied
    return [(l_alpha, l_beta) for l_alpha, l_beta in q if all(a + b == o  for o,a,b in zip(occ,l_alpha,l_beta)) ]

@timing
def sparce1(occ, n_alpha, n_beta):

    n_orbital = len(occ)

    occ_if = []
    occ_f = []
    for i, o in enumerate(occ):
        if o:
            occ_if.append(i)
            occ_f.append(o)

    #Add back the zeros
    def padd(l_occ):
        l = [0]*n_orbital
        for idx, o in zip(occ_if,l_occ):
            l[idx] = o
        return tuple(l)


    p =  list(product([0,1],repeat=len(occ_f)))
    p_alpha = filter(lambda o: sum(o) == n_alpha, p )
    p_beta =  filter(lambda o: sum(o) == n_beta, p )

    # All combinaison
    q = list(product(p_alpha, p_beta))

    # No filter corresponding to the correct occupation
    return [ tuple(map(padd,alpha_beta)) for alpha_beta in q if all( (a+b) == o  for o,a,b in zip(occ_f, *alpha_beta)) ] 

@timing
def double1(occ, n_alpha, n_beta):

    n_orbital = len(occ)

    occ_f = []
    occ_if = []

    for i,o in enumerate(occ):
        if o !=2:
            occ_f.append(o)
            occ_if.append(i)

    #Add back the zeros
    def padd(l_occ):
        l = [1]*n_orbital
        for idx, o in zip(occ_if,l_occ):
            l[idx] = o
        return tuple(l)

    n_double = len(occ) - len(occ_f)

    p =  list(product([0,1],repeat=len(occ_f)))
    # Filter alpha
    p_alpha = filter(lambda o: sum(o) == (n_alpha - n_double), p )
    p_beta =  filter(lambda o: sum(o) == (n_beta - n_double), p )

    # All combinaison
    q = product(p_alpha, p_beta)

    # No filter more than double occupied
    return [ tuple(map(padd,alpha_beta)) for alpha_beta in q if all( (a+b) == o  for o,a,b in zip(occ_f, *alpha_beta)) ]


@timing
def double1sparce1(occ, n_alpha, n_beta):

    n_orbital = len(occ)

    occ_if = []
    occ_f = []

    n_double = 0

    l_base = [None] * n_orbital

    for idx, o in enumerate(occ):
        if o == 1:
            occ_if.append(idx)
            occ_f.append(o)
        elif o == 0:
            l_base[idx] = 0
        elif o == 2:
            l_base[idx] = 1
            n_double += 1
             
    #Add back the zeros
    def padd(l_occ):

        l = l_base[:]
        for idx, o in zip(occ_if,l_occ):
            l[idx] = o

        return tuple(l)

    p =  list(product([0,1],repeat=len(occ_f)))
    p_alpha = filter(lambda o: sum(o) == (n_alpha - n_double), p )
    p_beta =  filter(lambda o: sum(o) == (n_beta - n_double), p )

    # All combinaison
    q = list(product(p_alpha, p_beta))

    # No filter corresponding to the correct occupation
    return [ tuple(map(padd,alpha_beta)) for alpha_beta in q if all( (a+b) == o  for o,a,b in zip(occ_f, *alpha_beta)) ]


@timing
def double1sparce1opt1(occ, n_alpha, n_beta =None):

    n_orbital = len(occ)

    occ_if = []
    n_double, n_orbital_f = 0, 0

    l_base = [None] * n_orbital

    for idx, o in enumerate(occ):
        if o == 1:
            occ_if.append(idx)
            n_orbital_f += 1
        elif o == 0:
            l_base[idx] = 0
        elif o == 2:
            l_base[idx] = 1
            n_double += 1

    c =  combinations(range(n_orbital_f), (n_alpha - n_double) )
   
    l_alpha = [0] * n_orbital_f

    l = []
    for sparce in c:
        dense = l_alpha[:]
    
        l_dense_base1 = l_base[:]
        l_dense_base2 = l_base[:]

        for i in sparce:
            dense[i] = 1
   
        for idx, o in zip(occ_if,dense):
            l_dense_base1[idx] = o
            l_dense_base2[idx] = not(o)

        l.append( (tuple(map(int,l_dense_base1)), tuple(map(int,l_dense_base2))))

    return l


from math import factorial

@timing
def double1sparce1opt1fortranlike1(occ, n_alpha, n_beta = None):

    n_orbital = len(occ)

    occ_if = [None]*n_orbital

    n_alpha_f, n_orbital_f = n_alpha, 0

    l_base = [None] * n_orbital

    for idx in range(n_orbital):
        o = occ[idx]
        if o == 1:
            occ_if[n_orbital_f] = idx
            n_orbital_f += 1
        elif o == 0:
            l_base[idx] = 0
        elif o == 2:
            l_base[idx] = 1
            n_alpha_f -= 1

    # Compute the number of det possible    
    d1 = 1
    for i in range(n_alpha_f):
         d1 *= (i+1)

    n1 = 1
    for i in range(n_orbital_f - n_alpha_f, n_orbital_f):
        n1 *= (i+1)
    
    n_det = n1 // d1 

    l_alpha = [0] * n_orbital_f
    l_det = [None] * n_det

    # Gen the permutation
    c =  list(combinations(range(n_orbital_f), (n_alpha_f) ))

    for idx_d in range(n_det):


        l_dense_base1 = l_base[0:n_orbital]
        l_dense_base2 = l_base[0:n_orbital]

        dense = [0] * n_orbital_f        
        sparce = c[idx_d]
        

        # n_alpha_f -> n_orbital_f
        for idx in range(n_alpha_f):
            i= sparce[idx]
            dense[i] = 1
            
        # n_orbtai_f -> n_orbital
        for i in range(n_orbital_f):
            o = dense[i]
            idx = occ_if[i]
            l_dense_base1[idx] = o
            l_dense_base2[idx] = not(o)

        l_det[idx_d] = ( tuple(l_dense_base1), tuple(l_dense_base2) )

    return l_det

import subprocess

from itertools import zip_longest
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def run_cf(t,  occ, n_alpha, n_beta):
    str_ = f"./{t}.out 1 {len(occ)} {n_alpha} {' '.join(map(str,occ))}"
    print(str_)
    s = subprocess.check_output(str_,shell=True, universal_newlines=True)
    l2 = []
    l = s.rstrip().split('\n')[2:]
    n_orbital = len(occ)

    n_int = n_orbital // 32 + 1
 

    for l in grouper(l,2*n_int):

        a = []
        for d in l[:n_int]:
            d = d.strip()
            d = '0' * (32 -len(d)) + d
            a.extend(d[::-1])

        b = []
        for d in l[n_int:]:
            d = d.strip()
            d = '0' * (32 - len(d)) + d
            b.extend(d[::-1])

        a = tuple([int(i) for i in a[:n_orbital]])
        b = tuple([int(i) for i in b[:n_orbital]])
        l2.append( (a,b) )
    return l2

@timing        
def c1(occ, n_alpha, n_beta): 
    return run_cf('c',occ, n_alpha, n_beta)

@timing
def f1(occ, n_alpha, n_beta):
        return run_cf('f',occ, n_alpha, n_beta)

def cmp_imp(f1: Callable, f2: Callable, n_up: int, n_down: int, seed: int):
    n_alpha, n_beta, occ = gen_occ(n_up,n_down, seed)
    str_ = f"occ:{occ} n_α: {n_alpha} n_β: {n_beta} seed: {seed}"
    print (str_)
    n = f1(occ, n_alpha, n_beta)
    o = f2(occ, n_alpha, n_beta)
    if sorted(n) != sorted(o):
        raise AssertionError(str_)
    print ()


def run_test(f1,f2,min_,max_):

    g = given(st.integers(min_value=min_, max_value=max_), st.integers(min_value=min_, max_value=max_), st.integers())

    def w(n_up,n_down, seed) -> Callable:
                cmp_imp(f1,f2,n_up,n_down, seed)

    g(w)()

if __name__ == '__main__':
    print ('Testing c')
    run_test(double1sparce1opt1, c1, 1,15)

    print ('Testing f')
    run_test(double1sparce1opt1, f1, 1,15)

