SWE/src/examples
================

Contains example programs that use SWE.

+ **swe_simple.cpp** A "simple" example that only runs on one core. Instead of the CPU it can also use the GPU for wave propagation.
+ **swe_mpi.cpp** Similar to the example above, but it can run on more the one node using MPI. If used with CUDA it requires one GPU per MPI task.
+ **swe_opengl.cpp** An example program that uses the OpenGL visualization.