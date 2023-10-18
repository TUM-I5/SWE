# Docker

## Build Docker Image
```shell
docker build -t swe .
```

## Run Docker image
```shell
docker run -it -v ${PWD}:/work --rm --privileged swe /bin/bash
```

A prebuilt Docker image is available in [Dockerhub](https://hub.docker.com/r/tumi5/swe).
