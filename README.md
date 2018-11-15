# Compile

```
>ifort -O2 occ2det.f90  -o f.out
>icc   -O2  occ2det.c  -o c.out
```


# Run
mode = 0, production / No output
mode = 1, testing / Output

```
>./$bin $mode $n_orbital $n_alpha $occ
```

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
