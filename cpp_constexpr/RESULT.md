gcc version 7.3.0 (Ubuntu 7.3.0-27ubuntu1~18.04) 

```bash
2018-11-18 23:20:54
Running ../../build/test/benchmark_test
Run on (8 X 3800 MHz CPU s)
CPU Caches:
  L1 Data 32K (x4)
  L1 Instruction 32K (x4)
  L2 Unified 256K (x4)
  L3 Unified 6144K (x1)
Load Average: 0.08, 0.10, 0.21
-------------------------------------------------------------
Benchmark                      Time           CPU Iterations
-------------------------------------------------------------
no_constexp_divmod      31653330 ns   31653214 ns         22
constexp_divmod         15038505 ns   15038445 ns         47
constexp_no_divmod      16026599 ns   16026545 ns         44
no_constexp_no_divmod   12506212 ns   12506245 ns         55
```