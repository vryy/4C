# BACI

BACI ("Bavarian Advanced Computational Initiative") is a parallel multiphsycis research code 
to address a plethora of physical problems by means of _computational mechanics_. 

Large parts of BACI are based on finite element methods (FEM), 
but alternative discretization methods such as discontinuous Galerkin methods (DG), 
particle methods and mesh-free methods have also been successfully integrated. 
The research software is implemented throughout in object-oriented programming (C++) 
using modern software design and is parallelized with MPI for distributed memory hardware architectures.

## Contents

1. [Getting Up and Running with BACI](#getting-up-and-running-with-baci)
   1. [First-Time Setup of Git](#first-time-setup-of-git)
   1. [Clone the Repository](#clone-the-repository)
   1. [Set Up the Environment](#set-up-the-environment)
   1. [Configure and Build](#configure-and-build)
   1. [Updating BACI](#updating-baci)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing](#contributing)
1. [License](#license)

## Getting Up and Running with BACI

### First-Time Setup of Git
You have to complete the following steps only once on your computer.
##### Your Identity
The first thing you should do in Git is setting your identity, i.e., your username and email address. This is important because every Git commit you create uses this information, and once a commit is finished this information is unchangeable.
Please set your username to your full name, i.e., first name followed by last name,
and your email address to your institute email address with the following commands:

```bash
git config --global user.name "Max Mustermann"
git config --global user.email mustermann@lnm.mw.tum.de
```
> **Note:** You may want to use a different name or email address for other projects your are working on. For that purpose, you can run the above commands without the `--global` option when you are in a project folder.

##### Your Text Editor
You can configure the default text editor that will be used whenever you need to write a message in Git.
The following command will set your default text editor to `kwrite`, a gui-based editor.

```bash
git config --global core.editor kwrite
```

> **Note:** Another popular choice is `vim`.

If you choose not to set a specific editor, Git will use your system’s default editor.

##### Check Your Settings
To confirm the correct setup of Git, you may check your configuration settings with:

```bash
git config --list
```
[↑ Contents](#contents)

### Clone the Repository

```bash
cd <someBaseDir>
mkdir <sourceDir>
git clone git@gitlab.lrz.de:baci/baci.git <sourceDir>
cd <sourceDir>
```

where `<someBaseDir>` is some directory on your machine and `<sourceDir>` will contain the BACI source code.
You can choose names and locations of these directories freely.

Your directory tree should look like the following:
```
<someBaseDir>/
  <sourceDir>
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
- [deal.II](www.dealii.org)

Often, a pre-compiled version of Trilinos and set of TPLs is available at your institute.
Look into the build configuration files in ```buildconfig/``` or ask your colleagues for further information.

[↑ Contents](#contents)

### Configure and Build

#### Create the Build Directory

BACI enforces an out-of-source build, i.e. your build directory may not be located inside the source code directory.

```bash
cd <someBaseDir>
mkdir <buildDir>
cd <buildDir>
```

where `<buildDir>` is your build directory.

#### Configure

Run

```bash
cd <someBaseDir>/<buildDir>
<someBaseDir>/<sourceDir>/do-configure | tee config$(date +%y%m%d%H%M%N).log
```

> **Note:**  When you see `command |& tee something$(date +%y%m%d%H%M%N).log`, that is just a means of running a command and sending the output both to the screen and to a timestamped log file.  This is by no means necessary, but if you run into problems, having these timestamped log files can be quite useful in debugging what's gone wrong.

#### Build

```bash
make -j <numProcs> |& tee make$(date +%y%m%d%H%M%N).log
```

where `<numProcs>` is the number of processors you want to use.

> **Note:**  After the first build, it is not always necessary to rerun the configure script &mdash; only the `make` command is required.  Reconfiguring is required when new files have been added and no changes are made to the `CMakeLists.txt` files.  If changes are made to a `CMakeLists.txt` file, then calling `make` will *automatically* reconfigure as part of the build process.

#### Run the Tests

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
cd <someBaseDir>/<sourceDir>
git checkout master
git pull
```

[↑ Contents](#contents)

## Where to Ask Questions

If you need help with BACI, feel free to ask questions by [creating a GitLab issue](https://gitlab.lrz.de/baci/baci/issues).  Use an issue template from the dropdown menu to pre-populate the *Description* field, giving you instructions on submitting the issue.

[↑ Contents](#contents)

## Contributing

If you're interested in contributing to BACI, we welcome your collaboration.  Please read [our contributing guidelines](https://gitlab.lrz.de/baci/baci/blob/master/CONTRIBUTING.md) for details on our workflow, submitting merge-requests, etc.

[↑ Contents](#contents)

## License

ADD A LICENSE!!

[↑ Contents](#contents)
