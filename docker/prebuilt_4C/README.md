>**Disclaimer:** this docker image is new. So, the usability might still be suboptimal.

This `Dockerfile` creates an docker image with a ready to use 4C executable.

You can just start the docker image with
```
docker run --interactive --tty ghcr.io/4c-multiphysics/4c:latest
```
or build it yourself with
```
cd docker/prebuilt_4C
docker build --tag 4c:latest .
```

Then you can run one of the input files inside the docker image with
```
mpirun -np 2 /home/user/4C/build/4C /home/user/4C/tests/input_files/<some input>.dat /home/user/output
```

If you want to use your own input file you can either use a volume mount (`--mount`) or copy the file into the docker container (`docker cp`).
