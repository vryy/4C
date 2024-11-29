# Docker

This folder contains Dockerfiles for different configurations. The `Dockerfile` in this folder is the main Dockerfile to
build the project with the current dependencies. Other configurations reside in different subfolders.

## Build docker images locally

Build the Docker image for the dependencies with the following command. **Note** that we need the correct build
context in order to copy the installation scripts. So, run the `docker build` in the projects root folder:

```bash
cd <project_root>
docker build --tag 4c-dependencies --file docker/Dockerfile .
```

## How to update the docker image for the Github Actions

1. Make your changes to the Dockerfile or dependencies
1. Compute the sha1 of the docker and dependencies folder

```bash
cd <project_root>
./docker/compute_dependencies_hash.sh
```

1. Update the `FOUR_C_DOCKER_DEPENDENCIES_HASH` variable in `docker.yml` and all mentions of the old hash (search and replace) in all workflows `.github/workflows` file (The hash should only occur in `container.image`).
1. Push a branch directly to the 4C repo (not your fork)
1. Open a PR (the automatically triggered workflows will fail)
1. Trigger the docker workflow manually on the branch in the 4C repo to build the docker image.
1. Update the branch such that the workflows are triggered again.
