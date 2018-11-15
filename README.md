# Compile

```
>ifort -O2 occ2det.f90  -o f.out
>icc   -O2  occ2det.c  -o c.out
```


# Run

```
>./$bin $mode $n_orbital $n_alpha $occ
```

mode = 0, production / No output
mode = 1, testing / Output


```
>time ./c.out 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1
n_int: 2
n_det: 346104

real    0m0.071s
user    0m0.051s
sys     0m0.009s

>time ./f.out 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1
n_int:           2
n_det:                346104

real    0m0.087s
user    0m0.075s
sys     0m0.011s
```

# Validation

Use ![hypothesis](https://github.com/HypothesisWorks/hypothesis) to generate random valide
occupation pattern and number of electrons.
We check that the C/Fortran implementation give the same result as a naive Python implementation.
```
def naive(occ, n_alpha, n_beta):

    n_orbital = len(occ)
    p =  list(product([0,1],repeat=n_orbital))
    # Filter the determinant with the correct number of alpha and beta
    p_alpha = filter(lambda o: sum(o) == n_alpha, p )
    p_beta =  filter(lambda o: sum(o) == n_beta, p )

    # All combinaison
    q = product(p_alpha, p_beta)

    # No filter the determinant who correspond to the required occupation
    return [(l_alpha, l_beta) for l_alpha, l_beta in q if all(a + b == o  for o,a,b in zip(occ,l_alpha,l_beta)) ]
```


```
./validation.py
```
