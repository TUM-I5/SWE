#pragma once

// Datatype for the type of data stored in the structures
#ifdef ENABLE_SINGLE_PRECISION
using RealType = float;
#define MY_MPI_FLOAT MPI_FLOAT
#else
using RealType = double;
#define MY_MPI_FLOAT MPI_DOUBLE
#endif
