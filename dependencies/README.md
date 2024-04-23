This folder contains the dependencies required to build the project. The current dependencies are in `current/`. To
install the dependencies just run

```
# Define the number of procs for building or use the default (4 procs)
export NPROCS=12
# Install each dependency
./current/<dependency>/install.sh <path/to/install/dir>
```

This folder also contains install scripts to build the Trilinos develop branch.
