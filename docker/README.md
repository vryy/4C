# Dockerfiles to build Docker images for Baci
Build the Docker image for the baci dependencies with the following command. **Note** that we need the correct build context in order to copy the installation scripts. So, run the `docker build` in the projects root folder
```
cd <project_root>
docker build --tag baci-dependencies --file docker/Dockerfile .
```
