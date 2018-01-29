[![Build Status](https://travis-ci.org/JBPennington/PLA.svg?branch=master)](https://travis-ci.org/JBPennington/PLA) [![codecov](https://codecov.io/gh/JBPennington/PLA/branch/master/graph/badge.svg)](https://codecov.io/gh/JBPennington/PLA)

# PLA

PLA (Portable Linear Algebra) is designed to be a portable, single header, MISRA compliant C, linear algebra library.

It's been built off of [datenwolf's linmath.h](https://github.com/datenwolf/linmath.h).

I cannot guarantee MISRA compliance at this time, but it is compliant to my knowledge.

## Cloning

Their is no need to clone, you could simply copy the header. If you'd like to clone and test the lib. Execute the following.

```
git clone --recursive https://github.com/JBPennington/PLA.git
```

-- or to clone as a submodule and integrate as a static lib --

```
git submodule add https://github.com/JBPennington/PLA.git
git submodule update --init --recursive
```

## Building

The only reason to build the repository is to test it.

To build, execute the following.

```
mkdir build; cd build;
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DTEST=ON; make;
```

A ctest executable "test_PLA" should've been made.

## Testing

At the moment, there should be 100% code coverage. To test, run the following...

```
ctest
```
-- or --
```
./bin/test_PLA
```
