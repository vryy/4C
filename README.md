# BACI

BACI ("Bavarian Advanced Computational Initiative") is a parallel multiphsycis research code 
to address a plethora of physical problems by means of _computational mechanics_. 

Large parts of BACI are based on finite element methods (FEM), 
but alternative discretization methods such as discontinuous Galerkin methods (DG), 
particle methods and mesh-free methods have also been successfully integrated. 
The research software is implemented throughout in object-oriented programming (```C++```) 
using modern software design and is parallelized with ```MPI``` for distributed memory hardware architectures.

## Contents

1. [Getting Up and Running with BACI](#getting-up-and-running-with-baci)
   1. [Clone the Repository](#clone-the-repository)
   1. [Set Up the Environment](#set-up-the-environment)
   1. [Set Up your Compiler Toolchain](#set-up-your-compiler-toolchain)
   1. [Configure and Build](#configure-and-build)
   1. [Updating BACI](#updating-baci)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing](#contributing)
1. [License](#license)

## Getting Up and Running with BACI

### Clone the Repositories

```bash
cd <someBaseDir>
git clone git@gitlab.lrz.de:baci/baci.git
cd baci
```

Your directory tree should look like the following:
```
<someBaseDir>/
  baci
```

[↑ Contents](#contents)

### Set Up the Environment

BACI heavily relies on the [Trilinos project](www.trilinos.org).

Some further third party libraries (TPLs) are mandatory, e.g.
- Parmetis
- SuiteSparse
- SuperLUDist
- Qhull

and some are optional, e.g.
- FFTW

Often, a pre-compiled version of Trilinos and set of TPLs is available at your institute.
Look into the build configuration files in ```buildconfig/``` or ask your colleagures for further information.

### Configure and Build

#### Create the Build Directory

BACI enforces an out-of-source build. 

```bash
cd ${WORKSPACE}
mkdir baci-${COMPILE_TYPE}-${BUILD_TYPE}
cd baci-${COMPILE_TYPE}-${BUILD_TYPE}
```

#### Configure

Run
```bash
cd ${WORKSPACE}/baci-${COMPILE_TYPE}-${BUILD_TYPE}
<someBaseDir>/baci/do-configure | tee config$(date +%y%m%d%H%M%N).log
```

> **Note:**  When you see `command |& tee something$(date +%y%m%d%H%M%N).log`, that is just a means of running a command and sending the output both to the screen and to a timestamped log file.  This is by no means necessary, but if you run into problems, having these timestamped log files can be quite useful in debugging what's gone wrong.

#### Build

```bash
make -j <numProcs> |& tee make$(date +%y%m%d%H%M%N).log
```

where `<numProcs>` is the number of processors you want to use.

> **Note:**  After the first build, it is not always necessary to rerun the configure script&mdash;only the `make` command is required.  Reconfiguring is required when new files have been added and no changes are made to the `CMakeLists.txt` files.  If changes are made to a `CMakeLists.txt` file, then calling `make` will *automatically* reconfigure as part of the build process.

#### Run the Unit Tests

To verify that the build was successful, run the minimal set of tests via
```bash
ctest -L minimal
```

or all tests via
```bash
ctest
```

[↑ Contents](#contents)

### Updating BACI

Any time you need to grab the latest from BACI:
```bash
cd ${WORKSPACE}/baci
git pull
```

[↑ Contents](#contents)

## Where to Ask Questions

If you need help with BACI, feel free to ask questions by [creating a GitLab issue](https://gitlab.lrz.de/baci/baci/issues).  Use an issue template to pre-populate the *Description* field, giving you instructions on submitting the issue.

[↑ Contents](#contents)

## Contributing

If you're interested in contributing to BACI, we welcome your collaboration.  Please read [our contributing guidelines](https://gitlab.lrz.de/baci/baci/blob/master/CONTRIBUTING.md) for details on our workflow, submitting merge-requests, etc.

[↑ Contents](#contents)

## License

ADD A LICENSE!!

[↑ Contents](#contents)
