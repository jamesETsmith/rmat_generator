# rmat_generator
A random matrix generator using recursive matrix (R-MAT) method from [Chakrabarti, Zhan, and Faloutsos](https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf).


## Build

```shell
cmake -B build
cmake --build build --parallel 6
ctest --test-dir build --parallel 6
```

## Benchmark

See the benchmark source files in `benchmark`. The built executables should all be in `<your_build_dir>/benchmark`. The command line interface for benchmarks should all support `--help` so run `<benchmark> --help` for up-to-date info about their command line interface. 