Stage0 += baseimage(image='ubuntu:22.04')
Stage0 += packages(ospackages=['python3','python3-pip',
    'build-essential', 'make', 'ninja-build',
    'autoconf','automake','libtool', 'git', 'less', 'wget', 'curl', 'bzip2', 'pkg-config', 'vim', 'unifdef',
    'texlive', 'texlive-full',
    'libopenmpi-dev', 'netcdf-bin', 'libnetcdf-dev', 'libnetcdf-c++4', 'libnetcdf-c++4-dev', 'libnetcdf-cxx-legacy-dev',
    'gdb', 'valgrind', 'cppcheck', 'ccache', 'graphviz',
    'clang-format', 'clang-tidy', 'doxygen', 'iwyu'])
Stage0 += gnu()
Stage0 += cmake(eula=True)
Stage0 += shell(commands=[
    'python3 -m pip install --upgrade pip matplotlib numpy jupyterlab hpccm',
])
Stage0 += environment(variables={
    'OMPI_ALLOW_RUN_AS_ROOT':'1',
    'OMPI_ALLOW_RUN_AS_ROOT_CONFIRM':'1',
    'DISPLAY':'host.docker.internal:0.0'
})
