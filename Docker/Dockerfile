FROM ubuntu:22.04
LABEL SWE=https://github.com/TUM-I5/SWE

ENV FORCE_UNSAFE_CONFIGURE=1 \
DEBIAN_FRONTEND=noninteractive \
TZ=UTC \
DISPLAY=host.docker.internal:0.0 \
OMPI_ALLOW_RUN_AS_ROOT=1 \
OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y --no-install-recommends --fix-missing \
software-properties-common \
apt-utils \
build-essential \
bash \
wget \
git \
curl \
ssh \
subversion \
tar \
unzip \
patch \
gzip \
bzip2 \
file \
gnupg \
coreutils \
mercurial \
nano \
pkg-config \
tree \
flex \
automake \
make \
cmake \
cmake-curses-gui \
ccache \
cppcheck \
libtool \
gcovr \
vim \
less \
meld \
gedit \
gdb \
strace \
linux-tools-common \
linux-tools-generic \
valgrind \
unifdef \
ninja-build \
meson \
graphviz \
iwyu \
doxygen \
xorg-dev \
mesa-utils \
mesa-common-dev \
freeglut3-dev \
ca-certificates \
libxml2-dev \
libbz2-dev \
libxslt1-dev \
libelf-dev \
libtbb-dev \
libtbb2 \
libsimde-dev \
libcurl4-openssl-dev \
libgl-dev \
libgl1-mesa-dev \
libglu1-mesa-dev \
libglvnd-core-dev \
libyaml-cpp-dev \
libnetcdf-c++4 \
libnetcdf-c++4-dev \
libnetcdf-cxx-legacy-dev \
libnetcdf-dev \
libopenmpi-dev \
netcdf-bin \
nodejs \
gcc \
g++ \
gfortran \
libgomp1 \
libnuma-dev \
libnuma1 \
g++-multilib \
gcc-multilib \
zlib1g \
zlib1g-dev \
python3 \
python3-dev \
python3-venv \
python3-pip && \
apt-get clean -y && \
apt-get autoremove -y && \
apt-get autoclean -y && \
rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade \
pip \
clang-format \
clang-tidy
