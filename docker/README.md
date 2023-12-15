# Baci in Docker
This folder contains Dockerfiles for different configurations. The `Dockerfile` in this folder is the main Dockerfile to build Baci with the current dependencies. Other configurations reside in different subfolders.

## Build docker images locally

Build the Docker image for the baci dependencies with the following command. **Note** that we need the correct build context in order to copy the installation scripts. So, run the `docker build` in the projects root folder:
```bash
cd <project_root>
docker build --tag baci-dependencies --file docker/Dockerfile .
```

## How to update the docker image

1. Make your changes to the Dockerfile or dependencies
1. Compute the sha1 of the docker and dependencies folder
```bash
cd <project_root>
find dependencies docker -not -wholename '*/trilinos_develop/*' -not -name 'README.md' -type f -exec sha1sum {} \; | sort | sha1sum | cut -c -8
```
1. Update `BACI_DOCKER_DEPENDENCIES_HASH` in the `.gitlab-ci.yml` file.
1. Run a manual pipeline with `BACI_DOCKER_BUILD_IMAGES: True`
