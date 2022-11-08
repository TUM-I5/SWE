include(FetchContent)

FetchContent_Declare(
    SWE-Solvers
    GIT_REPOSITORY https://github.com/TUM-I5/SWE-Solvers.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(SWE-Solvers)
