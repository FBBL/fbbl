




# FBBL - File-Based BKW for LWE 
[![Main Branch](https://img.shields.io/badge/-main:-gray.svg)](https://github.com/FBBL/fbbl/tree/main)&nbsp;&nbsp;&nbsp;&nbsp;[![Build Status](https://travis-ci.com/FBBL/fbbl.svg?token=xdMAmm6EEEu8xxxUqE6x&branch=master)](https://travis-ci.com/FBBL/fbbl)

**FBBL** (File-Based BKW for LWE) is a library that implements the [Blum-Kalai-Wasserman](https://arxiv.org/pdf/cs/0010022.pdf) (**BKW**) algorithm for solving relatively large LWE instances. It is part of the outcome of the research work presented with the paper [Making the BKW Algorithm Practical for LWE]() by Alessandro Budroni, Qian Guo, Thomas Johansson, Erik Mårtensson and Paul Stankovski Wagner.

The primary characteristic is to use a *file-based memory management* in order to avoid RAM constraints when handling large number of samples. It runs both on UNIX and Windows environments and it has no external dependencies except for [FFTW](http://fftw.org/) for solving Fast Fourier Transforms (FFT).

## Algorithm
### Input Problems and Sample Amplification

FBBL is able to generate LWE instances with round Gaussian error distribution and with an arbitrary number of samples. In case of a given external problem with limited number of samples (such as the case of the [TU Darmstadt LWE Challenges](https://www.latticechallenge.org/lwe_challenge/challenge.php)), **Sample Amplification** is implemented to generate the necessary number of samples. See the files in `/testVectors` for formatting your problem.

### Reduction Phase
FBBL supports most of the proposed BKW reduction steps for solving LWE problems, including:
- [plain-BKW](https://eprint.iacr.org/2012/636.pdf)
- [lazy modulus switching](https://eprint.iacr.org/2014/019.pdf) (LMS)
- [coded-BKW](https://eprint.iacr.org/2016/310.pdf)
- smooth-LMS

All the above support both **LF1** and **LF2** combining techniques.

**Coded-BKW** includes `(2,1), (3,1), (4,1)` and `(4,2)` codes. The latter is the concatenation of two `(2,1)` codes. The only modulos supported are `q = 101,631,1601,2053,16411.`

**Smooth-LMS** is the most general and advanced among them. It covers the funtionalities of plain-BKW and LMS, and it supports the **Unnatural Selection** tecnique. 

### Guessing Phase
For the guessing phase, FBBL supports both the **FFT**-based approach and the newly introduced (LINK) Fast-Walsh-Hadamard Transform (**FWHT**)-based approach. It is possible to install the library using only the latter and therefore with no external dependencies. The FFT solver can solve up to 3 positions at time while the FWHT solver has been tested to solve up to 33 (reduced) positions.

Both methods can run in combination with bruteforce to guess additional positions.

## Library
### Tests and Examples
The tests in `/test` are good examples to understand how the library works. They are all set to use `n=10, q=101` as parameters. The files in `/examples` are back-up files that will be converted into tests.

### Code formatting
The code is formatted using [Astyle](http://astyle.sourceforge.net/). Please run `format.sh` before making a PR or pushing your code.

## Authors
- [Alessandro Budroni](https://github.com/AlessandroBudroni)
- [Erik Mårtensson](https://github.com/ErikMaartensson)
- [Paul Stankovski](https://github.com/werekorren)


## How to Build on Linux and Mac

FBBL is a [Cmake](https://cmake.org/) project. Therefore one needs to install cmake to build the library.

### Dependencies

The only external library required is [FFTW](http://fftw.org/), used in the FFT-based guessing phase. It is possible to install FBBL witohut the support of the FFT solver and therefore using only the FWHT. In case one wants to have the FFT support, the library must be installed with the following three different settings:
1. `./configure --prefix=your_fftw_path && make && make install && make clean`
2. `./configure --prefix=your_fftw_path --enable-float && make && make install && make clean`
3. `./configure --prefix=your_fftw_path --enable-long-double && make && make install`

The location `YOUR_FFTW_PATH` will be used later when linking FFTW to FBBL. If not specified, the library will be installed in `/usr/local/`.

### Building Instructions

Customize the following `cmake` command to build the library:
```
mkdir target
cd target
cmake .. \
	-D PATH_PREFIX_A=your_path_for_data_A  \
	-D PATH_PREFIX_B=your_path_for_data_B  \
	-D FFTW_PREFIX=your_fftw_path \
	-D MAX_N=your_lwe_dimension \
	-D MAX_NUM_SAMPLES=your_max_num_samples \
	-D CMAKE_INSTALL_PREFIX=your_installation_path
make
```
Some explanation about the above `cmake` arguments:
- `PATH_PREFIX_A`: Location A for storing files, default: *${CMAKE_BINARY_DIR}/Testing*
- `PATH_PREFIX_B`: Location B for storing files, default: *${CMAKE_BINARY_DIR}/Testing*
- `MAX_N`: Max dimension n of the LWE problem, default: *10* (used by tests)
- `MAX_NUM_SAMPLES`: Max number of samples, default: *1000000*  (used by tests)
- `FFTW_PREFIX`: fftw installation path, default: `/usr/local`
- `CMAKE_INSTALL_PREFIX`: Install directory used by `make install`. default `/usr/local/`


Other `cmake` arguments:
- `EARLY_ABORT_LOAD_LIMIT_PERCENTAGE`: Early abort load limit percentage, default: *96*
- `MIN_STORAGE_WRITER_CACHE_LOAD_BEFORE_FLUSH`: Minimum storage writer cache load before flush, default: *25*
- `SAMPLE_DEPENDENCY_SMEARING`: Keep the number of samples to be somehow constant through the steps, default: *1*
- `CMAKE_BUILD_TYPE`: Compiler flags mode, default: *Release*, other possibilities: *Test* and *Coverage*
- `CMAKE_VERBOSE_MAKEFILE`: Cmake verbose output, default: *FALSE*
- `BUILD_TESTING`: Build tests, default: *ON*
- `BUILD_FFT`: Build FFT support (fftw library required), default: *ON*

The command 
```
make test
```
run some instances of the algorithm with different settings. The samples generated and the auxiliary data are stored in `PATH_PREFIX_A`. The following command remove such data.
```
make dataclean
```

### Coverage
Coverage of the code is produced using [lcov](http://ltp.sourceforge.net/coverage/lcov.php) with the command
```
make coverage
```
available only when `CMAKE_BUILD_TYPE=Coverage`.

### Install

Install the library in `your_installation_path` with the following command.
```
sudo make install
``` 

## HOW TO USE
After installation, compile (for example using GCC) your program with
```
gcc example.c -I your_installation_path/include -L your_installation_path/lib -fbbl -o example
```

## TODO/Wish List
- add build option and guidelines for Windows
- implement [coded-BKW with Sieving](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_12) reduction step
- implement a new file-based sample handling memory technique
- more tests
- code documentation

Pull requests are always welcome!