## Description and Context

This merge request updates the supported Trilinos version to X.X.X (git hash `abcdefg123456`).

## Checklist
- [ ] Copy the Trilinos develop scripts under `dependencies/trilinos_develop` to the current folder `dependencies/current/`. Don't forget to explicitly add the Trilinos hash in the Trilinos install script.
- [ ] Update the sha1 of the docker image and start a manual pipeline to build the docker image, see [docker/README.md](../../docker/README.md). Don't forget to set the new hash in the test configuration file and add it to the test configuration file (`tests/testconfig/.gitlab-ci.yml`).
- [ ] Update the Trilinos installation for LNM workstations
- [ ] Update the Trilinos installation for LNM clusters
- [ ] Update the Trilinos installation for IMCS workstations
- [ ] Update the Trilinos installation for IMCS clusters
- [ ] Update the Trilinos path in the IMCS preset files
