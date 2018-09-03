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
   1. [Set-up Git](#set-up-git)
   1. [Set-up LRZ GitLab](#set-up-lrz-gitlab)
   1. [Clone the Repository](#clone-the-repository)
   1. [Set Up the Environment](#set-up-the-environment)
   1. [Configure and Build](#configure-and-build)
   1. [Updating BACI](#updating-baci)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing](#contributing)
1. [License](#license)

## Getting Up and Running with BACI

### Set-up Git

1. Set your username to your full name, i.e., first name followed by last name,
and your email address to your institute email address with the following commands:

    ```bash
    git config --global user.name "<Firstname> <Lastname>"
    git config --global user.email <instituteEmailAddress>
    ```

1. Set a default text editor that will be used whenever you need to write a message in Git. To set `kwrite` as your default text editor, type: 

    ```bash
    git config --global core.editor kwrite
    ```

    > **Note:** Another popular choice is `vim`.

Our Wiki provides a [detailed setup guide for your local git configuration](https://gitlab.lrz.de/baci/baci/wikis/Set-up-Git).

[↑ Contents](#contents)

### Set-up LRZ GitLab

1. Register an account on [LRZ GitLab](www.gitlab.lrz.de).  
   **Important:** Choose a recognizable user name. It is recommended to set it to: first letter of first name followed by last name, all lowercase, e.g., Max Mustermann -> mmustermann.

    > **Note:** Your username is a unique namespace related to your user ID. Changing it can have unintended side effects. See [how redirects will behave](https://gitlab.lrz.de/help/user/project/index.md#redirects-when-changing-repository-paths) for details.

1. Go to your GitLab profile settings and update your profile settings, in particular your
    * First and last name
    * Institute Email Address
1. Select proper notification settings. We recommend *Watching* or *On mention* to guarantee that you don't miss any important developments and discussions.
1. Add your public SSH key found in `~/.ssh/id_rsa.pub` to your user profile.

Our Wiki provides [detailed setup instructions for your GitLab account](https://gitlab.lrz.de/baci/baci/wikis/Set-up-and-Configure-your-GitLab-Account).

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
<someBaseDir>/<sourceDir>/do-configure --config=<path/to/build-configuration-file.config> | tee config$(date +%y%m%d%H%M%N).log
```

> **Note:**  When you see `command |& tee something$(date +%y%m%d%H%M%N).log`, that is just a means of running a command and sending the output both to the screen and to a timestamped log file.  This is by no means necessary, but if you run into problems, having these timestamped log files can be quite useful in debugging what's gone wrong.

A build configuration file needs to be passed to the configure script via the command line argument `--config`, as indicated above.
Configuration files for a bunch of supported system environments are located in `<someBaseDir>/<sourceDir>/buildconfig/`.

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

If you're interested in contributing to BACI, we welcome your collaboration. Before your start, configure your [local Git](#set-up-git) and your [LRZ GitLab account](#set-up-lrz-gitlab). Read [our contributing guidelines](https://gitlab.lrz.de/baci/baci/blob/master/CONTRIBUTING.md) carefully for details on our workflow, submitting merge-requests, etc.

[↑ Contents](#contents)

## License

ADD A LICENSE!!

[↑ Contents](#contents)
