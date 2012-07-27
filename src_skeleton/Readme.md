Teaching skeleton
===

# CUDA skeleton
This directory holds the CUDA skeleton for the *Gene Golub SIAM Summer School 2012*.

There are two CUDA kernels in the file SWE_WavePropagationBlockCuda_kenels.cu, which need to be implemented:
1. void computeNetUpdatesKernel([...])
2. void updateUnknownsKernel([...])

A C++-reference is available in the file [SWE_WavePropagationBlock.cpp](https://github.com/TUM-I5/SWE/blob/master/src/SWE_WavePropagationBlock.cpp), remark: The sizes of the arrays differs from the planned CUDA implementation. This is outlined in more deatail within the constructor of [SWE_WavePropagationBlockCuda][https://github.com/TUM-I5/SWE/blob/3f9a316d196005d39496ce7231a57c6cf3961ec3/src/SWE_WavePropagationBlockCuda.cu#L52).

The [CUDA version of the f-wave solver](https://github.com/TUM-I5/swe_solvers/blob/master/src/solver/FWaveCuda.h) should be used for the computation of the net-updates.

Theres a tested CUDA implementation of the kernels implemented as well in the src directory of SWE [src/WE_WavePropagationBlockCuda_kernels.cu](https://github.com/TUM-I5/SWE/blob/master/src/SWE_WavePropagationBlockCuda_kernels.cu). Nevertheless we recommend to work with the C++-implementation as reference only due to the improved learning effect.

# Main file
An example main file is located at [src/examples/swe_wavepropagation.cpp](https://github.com/TUM-I5/SWE/blob/master/src/examples/swe_wavepropagation.cpp). Details about the compile and linking process can be found in the corresponding [SCons-script](https://github.com/TUM-I5/SWE/blob/master/src/SConscript), which selects the source files.
