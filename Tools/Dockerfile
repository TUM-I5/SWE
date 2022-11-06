FROM ubuntu:22.04

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        bzip2 \
        ccache \
        clang-format \
        clang-tidy \
        cppcheck \
        curl \
        doxygen \
        gdb \
        git \
        graphviz \
        iwyu \
        less \
        libopenmpi-dev \
        libtool \
        make \
        ninja-build \
        pkg-config \
        python3 \
        python3-pip \
        unifdef \
        valgrind \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran && \
    rm -rf /var/lib/apt/lists/*

# CMake version 3.22.2
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        make \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2-linux-x86_64.sh && \
    mkdir -p /usr/local && \
    /bin/sh /var/tmp/cmake-3.22.2-linux-x86_64.sh --prefix=/usr/local --skip-license && \
    rm -rf /var/tmp/cmake-3.22.2-linux-x86_64.sh
ENV PATH=/usr/local/bin:$PATH

RUN python3 -m pip install --upgrade pip matplotlib numpy jupyterlab hpccm

ENV DISPLAY=host.docker.internal:0.0 \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1


