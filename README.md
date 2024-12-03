<h1 align="center">
  4C
</h1>

<div align="center">

[![website](./utilities/assets/badges/website_badge.svg)](https://4C-multiphysics.org)
[![docs/rtd](./utilities/assets/badges/documentation_readthedocs.svg)](https://4c-multiphysics.github.io/4C/readthedocs/)
[![docs/doxygen](./utilities/assets/badges/documentation_doxygen.svg)](https://4c-multiphysics.github.io/4C/doxygen/)

</div>

<div align="center">

[![workflows/checkcode](https://github.com/4C-multiphysics/4C/actions/workflows/checkcode.yml/badge.svg?branch=main)](https://github.com/4C-multiphysics/4C/actions/workflows/checkcode.yml?query=branch%3Amain)
[![workflows/buildtest](https://github.com/4C-multiphysics/4C/actions/workflows/buildtest.yml/badge.svg?branch=main)](https://github.com/4C-multiphysics/4C/actions/workflows/buildtest.yml?query=branch%3Amain)
[![workflows/nightly_tests](https://github.com/4C-multiphysics/4C/actions/workflows/nightly_tests.yml/badge.svg?branch=main)](https://github.com/4C-multiphysics/4C/actions/workflows/nightly_tests.yml?query=branch%3Amain)
[![workflows/documentation](https://github.com/4C-multiphysics/4C/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/4C-multiphysics/4C/actions/workflows/documentation.yml?query=branch%3Amain)

</div>

> [!WARNING]
> We are currently migrating the project from GitLab. Links might be out-dated.

4C ("Comprehensive Computational Community Code") is a parallel multiphysics research code
to address a plethora of physical problems by means of _computational mechanics_.

Large parts of 4C are based on finite element methods (FEM),
but alternative discretization methods such as discontinuous Galerkin methods (DG),
particle methods and mesh-free methods have also been successfully integrated.
The research software is implemented throughout in object-oriented programming (C++)
using modern software design and is parallelized with MPI for distributed memory hardware architectures.

**Disclaimer**: 4C is developed for research purposes in the field of numerical method development.
It is not intended for any use beyond this purpose, and generally should not be used for any form of
safety-relevant or safety-critical calculations,
or for an application in association with physical products in particular.

## Required development tools

- **CMake:** 4C supports CMake configuration through CMake presets.
- **Ninja:** We strongly encourage to work with `ninja` (instead of `make`) for faster compile times.

## External dependencies

4C heavily relies on the [Trilinos project](https://trilinos.github.io).

Some further third-party libraries (TPLs) are mandatory, e.g.

- CMake (minimum version: 3.30)
- ParMETIS (recommended version: 4.0.3)
- SuiteSparse (recommended version: 5.4.0)
- SuperLUDist (recommended version: 7.2)
- Qhull (recommended version: 2012.1)
- CLN (recommended version: 1.3.4)

and some are optional, e.g.

- FFTW
- [ArborX](https://github.com/arborx/ArborX)
- [MIRCO](https://github.com/imcs-compsim/MIRCO/)

Maybe, a pre-compiled version of Trilinos and set of TPLs is available at your institute.
Look into the CMake presets in `presets/` or ask your colleagues for further information.

Some helper scripts to install these TPLs can be found in `dependencies/`.

Additional information can be found in
the [user documentation](https://4c-multiphysics.github.io/4C/readthedocs/4Csetup.html#external-dependencies).

## Getting up and running

### Create python virtual environment (required for code development in 4C)

For testing and active development, you need to create a python virtual environment once.
In the source directory,
execute:

```
./utilities/set_up_dev_env.sh
```

### Configure and Build

4C enforces an out-of-source build, i.e., the build directory may not be the top-level source directory.
Create a build directory `<path/to/build/dir>` of your choosing.

Navigate into the build directory and run

```bash
cmake --preset=<name-of-preset> <path/to/source/dir>
```

A preset name needs to be passed to cmake via the command line argument `--preset`.
Use `cmake <path/to/source/dir> --list-presets` to get a list of all available presets.
Define your own presets in a `CMakeUserPresets.json` file.

More information about CMake presets can be found in the [CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

#### Build

To build the code on `<numProcs>` processes, run

```bash
ninja -j <numProcs> full
```

#### Run the Tests

To verify that the build was successful, run the minimal set of tests via

```bash
ctest -L minimal
```

or all tests via

```bash
ctest
```

You can use the option `-j <num_threads>` to specify the number of threads to be used for parallel execution.

### Prepare and Run Simulations

After successfully building 4C, the executable `4C` is located in your build directory `<path/to/build/dir>/`.
It needs to be invoked together with an input (`.dat`) file via

```bash
<path/to/build/dir>/4C <jobName>.dat <outputName>
```

where `<jobName>` is the name of the simulation file and `<outputName>` denotes the name of the corresponding output
file(s).
A collection of working `.dat` files is located under `tests/input_files/` in the source code repository.

Input files (`*.dat`) can be generated through various mechanisms.
Please consult our [user guide](https://4c-multiphysics.github.io/4C/readthedocs/index.html) for further information and detailed instructions.

4C can write its simulation output in different formats.
The [user guide](https://4c-multiphysics.github.io/4C/readthedocs/index.html) outlines the different output format and their necessary steps to actually view the results.

## Where to Ask Questions

If you need help with 4C, feel free to ask questions
by [creating a GitHub issue](https://github.com/4C-multiphysics/4C/issues). Use an issue template to pre-populate the *Description* field, giving you instructions on submitting the issue.

For more general questions, you can reach us on the [4C Slack workspace](https://join.slack.com/t/4c-multiphysics/shared_invite/zt-1oi61jgdd-5tZuHku3Tb_BH5UBgojbpQ).

## Contributing

If you're interested in contributing to 4C, we welcome your collaboration.
Read [our contributing guidelines](https://github.com/4C-multiphysics/4C/blob/main/CONTRIBUTING.md) carefully for details on
our workflow, submitting pull requests, etc.

## Code of Conduct

All people and activities in and around 4C are subject to our [Code of Conduct](https://github.com/4C-multiphysics/4C/blob/main/CODE_OF_CONDUCT.md).

## How to cite 4C

Please cite 4C as follows:

```
4C: A Comprehensive Multi-Physics Simulation Framework, https://www.4c-multiphysics.org
```

You could use the following BibTeX entry:

```bibtex
@misc{4C,
  author       = {{4C}},
  title        = {{4C}: A {C}omprehensive {M}ulti-{P}hysics {S}imulation {F}ramework},
  howpublished = {\url{https://www.4c-multiphysics.org}},
  year         = {YEAR},
  note         = {Accessed: DATE}
}
```

We kindly ask you to also give credit to the individual methods and algorithms used in 4C.
References to the relevant publications can be found on the [4C website](https://4c-multiphysics.org) or throughout the source code.
If you need any assistance with finding suitable references,
please feel free to reach out in the [4C Slack workspace](https://join.slack.com/t/4c-multiphysics/shared_invite/zt-1oi61jgdd-5tZuHku3Tb_BH5UBgojbpQ).
