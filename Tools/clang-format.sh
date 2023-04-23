#!/bin/bash
for i in ../Source/**/*.c* ../Source/**/*.h* ../Source/*.c* ../Source/*.h* ../Tests/*.c* ../Tests/*.h* ; do echo $i && clang-format -i $i ;done
