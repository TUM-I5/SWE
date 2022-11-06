include(FetchContent)

FetchContent_Declare(
    swe-solvers
    GIT_REPOSITORY https://github.com/TUM-I5/SWE-Solvers.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(swe-solvers)
