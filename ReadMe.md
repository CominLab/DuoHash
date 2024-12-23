# DuoHash: fast hashing of spaced seeds with application to spaced k-mers counting



## Summary

Alignment-free genomic sequence analysis has facilitated high-throughput processing within numerous bioinformatics workflows. A central task in alignment-free applications is hashing $k$-mers, commonly used for indexing, querying, and fast similarity searches. Recently, spaced seeds—a specialized pattern designed to accommodate errors or mutations—have increasingly replaced $k$-mers, enhancing sensitivity in various applications. However, spaced seed hashing is computationally intensive, introducing significant delays.
This paper addresses the challenge of efficient spaced seed hashing and presents DuoHash,  a framework that enables the efficient computation of several hash functions. Our experimental results demonstrate that the proposed method substantially outperforms existing algorithms, achieving speedups of up to 11x. To illustrate practical utility, we further applied DuoHash to the problem of spaced $k$-mers counting.

## Methods
The **DuoHash** library provides two classes: **DuoHash** and **DuoHash_multi** for handling one or multiple spaced seeds, respectively. The methods of the first class are
- `GetEncoding_naive()`,
- `GetEncoding_FSH()`,
- and `GetEncoding_ISSH()`.

The methods of the second class are
- `GetEncoding_naive()`,
- `GetEncoding_FSH()`,
- `GetEncoding_ISSH()`,
- `GetEncoding_FSH_multi()`,
- `GetEncoding_MISSH_v1()`,
- `GetEncoding_MISSH_col()`,
- `GetEncoding_MISSH_col_parallel()`,
- and `GetEncoding_MISSH_row()`.

Both classes share the `PrintFASTA()` method for saving the resulting spaced k-mers to a file and other methods for handling the various parameters.

Each `GetEncoding_<...>()` method has four implementations. The first is for the extraction of spaced k-mer and their encoding only, the second allows post-processing of encodings to calculate forward and reverse hashing, the third allows post-processing of encodings for conversion into strings, and the fourth combines the two previous options.


## Installation
Make sure CMake is installed on the system.

Download the repository using
```shell
$ git clone https://github.com/leonardoGemin/DuoHash.git
```
and build the library with 
```shell
$ make build
```

This will install `build/libDuoHash.a` in the project's directory.


## Usage
To use **DuoHash** in a C++ project:
- Import DuoHash in the code using `#include <DuoHash.h>`
- Add the `include` directory (pass `-I./include` to the compiler)
- Link the code with `libDuoHash.a` (pass `-L./build -lDuoHash` to the compiler)
- Compile your code with `g++-13`, `-std=c++0x` (and preferably `-O3`), and `-fopenmp` enabled


## Example
Compile `example/main.cpp` file with
```shell
$ cd example
$ g++-13 -std=c++0x -O3 -fopenmp -I../include -L../build -lDuoHash -o main main.cpp
```


## Citation

DuoHash: fast hashing of spaced seeds with application to spaced k-mers counting

Leonardo Gemin, Cinzia Pizzi, and Matteo Comin

Accepted at ICCABS 2025


## Thesis
Link to Gemin Master Thesis: [Gemin_Leonardo.pdf](https://thesis.unipd.it/retrieve/99e7ee7c-1348-45f6-8a02-467bae6b0dbc/Gemin_Leonardo.pdf)
