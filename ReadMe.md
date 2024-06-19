# DuoHash: Improving Spaced k-mer Extraction and Hash Encoding for Bioinformatics Applications


## Methods
The **DuoHash** library provides two classes: **DuoHash** and **DuoHash_multi** for handling one or multiple spaced seeds, respectively. The methods of the first class are `GetEncoding_naive()`, `GetEncoding_FSH()`, and `GetEncoding_ISSH()`. The methods of the second class are `GetEncoding_naive()`, `GetEncoding_FSH()`, `GetEncoding_ISSH()`, `GetEncoding_FSH_multi()`, `GetEncoding_MISSH_v1()`, `GetEncoding_MISSH_col()`, `GetEncoding_MISSH_col_parallel()`, and `GetEncoding_MISSH_row()`. Both classes share the `PrintFASTA()` method for saving the resulting spaced k-mers to a file and other methods for handling the various parameters. Each `GetEncoding_<...>()` method has two implementations. The first is for the method itself, and the second allows for post-processing of the encoding through three possible functions (`GetHashes`, `GetSpacedKmer`, and `GetBoth`).


## Installation
Make sure CMake is installed on the system.

Download the repository using `git clone https://github.com/[...]` and build the library with 
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
