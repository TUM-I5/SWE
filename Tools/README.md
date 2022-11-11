# Dockerfile

## Generate Dockerfile using HPC Container Maker (HPCCM)
You can use (and adapt) `docker_recipe.py` to generate a Dockerfile. At first, you need to install `hpccm`:

```shell
python3 -m pip install --upgrade hpccm --user
```

Then, you can execute following shell script:

```shell
hpccm --recipe docker_recipe.py --format docker > Dockerfile
```
A basic Dockerfile has already been generated.

## Build Docker image
As a next step, you need to build the Docker image using the previously generated Dockerfile:

```shell
docker build -t swe .
```

This might take a while, since the full LaTeX distribution is also installed.

## Run Docker image
In order to run the Docker image as a container you need to execute the following command:

```shell
docker run \
  -it \
  --rm \
  --gpus=all \
  --privileged \
  --cap-add=SYS_PTRACE \
  --security-opt label=disable \
  --security-opt seccomp=unconfined \
  --security-opt apparmor=unconfined \
  --ipc=host \
  --network=host \
  --ulimit memlock=-1 \
  --ulimit stack=67108864 \
  -u $(id -u):$(id -g) \
  -v $(pwd):/work \
  -w /work \
  swe \
  /bin/bash
```

If you are using Windows, then you need to execute this command instead:

```shell
docker run \
  -it \
  --rm \
  --gpus=all \
  --privileged \
  --cap-add=SYS_PTRACE \
  --security-opt label=disable \
  --security-opt seccomp=unconfined \
  --security-opt apparmor=unconfined \
  --ipc=host \
  --network=host \
  --ulimit memlock=-1 \
  --ulimit stack=67108864 \
  -v <your\work\directory>\:/work \
  swe \
  /bin/bash
```
