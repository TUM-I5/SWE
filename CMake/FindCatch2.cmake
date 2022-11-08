include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.1.1
)
FetchContent_MakeAvailable(Catch2)

target_disable_static_analysis(Catch2)
target_disable_static_analysis(Catch2WithMain)
