# TODO - Group 1 (Tasks, Results & DOCS)

This is a list of what Group 1 did as tasks for intermediate documentation purposes:
1. [x] Move `BoundaryEdge.hpp` to Types folder - Types cleaning
2. [ ] Reduction of warnings `Wuseless-cast`, `-Wsign-conversion`
   1. [ ] After the fixed propose a merge request to github original repo
3. Nsight Compute Guided Analysis
4. Check solver Sparsity
5. Multi stream execution IMPORTANT
6. Unified Memory
7. what did you do in `.zshrc` to make it work?
8. how to compile and run the code in side the gpu cluste rlgon node-> x-login0 -> singularity shell...
9. UMA results:
10.   - 86 = sm_86 = Ampere 2nd gen (RTX 30xx)
  - `~/tooling/cuda-samples/Samples/6_Performance/UnifiedMemoryPerf`
  - CUDA Systems Integration, Unified Memory, CUDA Streams and Events, Pinned System Paged Memory:
  ```
  GPU Device 0: "Ampere" with compute capability 8.6

  Running ........................................................

  Overall Time For matrixMultiplyPerf 

  "UMhint",   // Managed Memory With Hints
  "UMhntAs",  // Managed Memory With_Hints Async
  "UMeasy",   // Managed_Memory with No Hints
  "0Copy",    // Zero Copy
  "MemCopy",  // USE HOST PAGEABLE AND DEVICE_MEMORY
  "CpAsync",  // USE HOST PAGEABLE AND DEVICE_MEMORY ASYNC
  "CpHpglk",  // USE HOST PAGELOCKED AND DEVICE MEMORY
  "CpPglAs"   // USE HOST PAGELOCKED AND DEVICE MEMORY ASYNC

  Printing Average of 20 measurements in (ms)
  Size_KB  UMhint UMhntAs  UMeasy   0Copy MemCopy CpAsync CpHpglk CpPglAs
  4         0.213   0.222   0.333   0.016   0.032   0.026   0.032   0.026
  16        0.234   0.253   0.467   0.028   0.043   0.046   0.054   0.044
  64        0.323   0.357   0.834   0.113   0.119   0.093   0.088   0.078
  256       0.566   0.592   1.163   0.469   0.282   0.267   0.248   0.235
  1024      2.516   1.967   2.540   2.595   1.134   1.096   0.976   0.961
  4096      6.533   5.672   8.640  13.487   4.846   4.806   4.298   4.283
  16384    29.209  26.390  38.098  93.926  22.238  22.168  21.201  21.265
  ```
  - `CUDA_MALLOC`
[x]  16. mention how the validation detected a bug in the SIMD XXX with diff scenarios
[x]  17. mention the optimization error due to -ffast-math -fassociative-math -freciprocal-math -fno-signed-zeros, as it works fine in DEBUG mode [](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html#floating-point-math-is-not-associative) [](https://developer.download.nvidia.com/assets/cuda/files/NVIDIA-CUDA-Floating-Point.pdf?t=eyJscyI6InJlZiIsImxzZCI6IlJFRi1kb2NzLm52aWRpYS5jb20vIn0=)
18. mention how you validated using paraview `plot over line` for all metrics h hu,,..
19. mention how we sat up gitlabrunner with intel compiler using docker and runnign on local server sam-legion
20. mention `CUDA_MEM_CPY` and how to use ccmake to enable and disable UMA and similar macros
21. Document in intrinsics the manual striding we did for dead cells in Float2D
22. Test inlining
23. clean cmake between differnt settings ENABLE_CUDA_UMA_ALLOCATOR
24. why pinned page slow?
25. streams update
26. Profiling:
   - [x] Intel - Advisor:
     - [x] For original SWE
       - [x] Roofline model
         - [x] Compute bound
         - [x] Plot - html
     - [ ] For vectorised SWE
         - [ ] Roofline model
           - [ ] Compute bound
           - [ ] Plot - html
     - [ ] Start tuning to reach 60-70% of the theoretical peak
       - [x] for now consider 41.6 GFLOPs
       - [ ] parallelize rows VS cols VS kernels(because of the blocking concept) to tune
       - [ ] Support OMP for our SIMD code
     - [x] 4 doubles  -> is it the size of the cache -> check the assignment 1 L1 cache ! @Nandini
       - ```bash
            t1221am@i22r07c05s03:~/project/SWE_sam/build>  getconf -a | grep CACHE
            LEVEL1_ICACHE_SIZE                 32768
            LEVEL1_ICACHE_ASSOC                8
            LEVEL1_ICACHE_LINESIZE             64
            LEVEL1_DCACHE_SIZE                 32768
            LEVEL1_DCACHE_ASSOC                8
            LEVEL1_DCACHE_LINESIZE             64
            LEVEL2_CACHE_SIZE                  262144
            LEVEL2_CACHE_ASSOC                 8
            LEVEL2_CACHE_LINESIZE              64
            LEVEL3_CACHE_SIZE                  18350080
            LEVEL3_CACHE_ASSOC                 20
            LEVEL3_CACHE_LINESIZE              64
            LEVEL4_CACHE_SIZE                  0
            LEVEL4_CACHE_ASSOC                 0
            LEVEL4_CACHE_LINESIZE              0
         ```
     - [x] Stream benchmark
        - [x] collect original bandwidth for cluster
     - [x] Calculate BW for the vectorized code + compare it to autovectorized code
       - [ ] Using sheet 3, we need to measure the bus bandwidth for CM2
         - [x] then consider the BW limit in our speedups
       - [ ] Compare the speedup and results of our with the literature(paper)
         - [ ] Theoratical peak
         - [ ] the hardware
         - [ ] GFLOPS calculation
         - [ ] Roofline model
     - [ ] LIKWID : (type of events : cache misses for L1, L2, L3 and total, bandwidth)
          - [ ] run perf counter for original SWE @Nandini
          - [ ] run perf counter for original omp auto vectorised SWE @Nandini
          - [ ] run perf counter for original vectorised SWE
          - [ ] run perf counter for original SWE with CUDA
          - [x] find what we did for Likwiid takd Q4 (assignment paper)
   - [ ]<***URGENT***> How to present our results/measurements for vectorization
     - [ ] In the paper, they reported "Million Riemann solves per second per core" (short: MRim/s per    core) for various matrix sizes.(Performance of fWave Solver for vectorized and non-vectorized). "Million element updates per second per core" (short: MElup/s) - each element update requires two Riemann solves and the respective numbers include the computational effort to update the unknowns based on the accumulated net updates.(Performance of entire SWE using vectorized and non-vectorized fWave solver)
     - [ ] finalize what we need to do for our code
   - [ ] Nvidia Insight Compute: @Sam
     - [ ] Kill unnecessary copies between device-host in UMA mode
     - [ ] Compare rows / columns / kernels of loading to GPU
     - [ ] GPU occupancy calculater (shows density)
     - [ ] Mem allocation type comparison between UMA and device mem
     - [ ] Profile + Tune using Intel adv
     - [ ] Plot warps and data in mem
     - [ ] data streams to optimize
27. XX URGENT XX:
    - [x] Fix the exit of vectorised code
        - [ ] Generate a new file from vectorised and run in Paraview -> fix werong resdiual first
    - [x] Disable both VTK and NetCDF
    - [ ] Support alligned memory
28. Document:
   - [x] Structure outline for the report(i.e. titles, neceassry sections, etc)
   - [ ] fetch all comments i added in the code as // FIXME and // NOTE and // TODO
   - [ ] add README and how to run the code for
     - [ ] vectorization
     - [ ] CUDA
     - [ ] any other options added by us
29. Unit and Regression Testing:
     - [ ] reestablish gitlab runner on Sams laptop
     - [ ] rerun the pipelines
       - [ ] vectorized code
       - [ ] cuda
     - [ ] Add more unit tests for verification + validations
       - [ ] Document

# References List
- https://developer.nvidia.com/blog/easy-introduction-cuda-c-and-c/
- https://developer.nvidia.com/blog/unified-memory-cuda-beginners/
- https://developer.nvidia.com/blog/how-implement-performance-metrics-cuda-cc/