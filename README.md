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
   1. [Clone the Repository](#clone-the-repository)
   1. [Set Up the Environment](#set-up-the-environment)
   1. [Configure and Build](#configure-and-build)
   1. [Updating BACI](#updating-baci)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing](#contributing)
   1. [First-Time Setup](#first-time-setup)
      1. [Set Up Git](#set-up-git)
      1. [Set Up GitLab](#set-up-gitlab)
1. [License](#license)

## Getting Up and Running with BACI

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

If you're interested in contributing to BACI, we welcome your collaboration. Before your start, please configure your [local Git](#set-up-git) and your [GitLab account](#set-up-lrz-gitlab) according to the instructions below. Next, please also read [our contributing guidelines](https://gitlab.lrz.de/baci/baci/blob/master/CONTRIBUTING.md) for details on our workflow, submitting merge-requests, etc.

[↑ Contents](#contents)

### First-Time Setup
Before your start contributing to BACI, please complete the following configuration steps. You have do to this only once.

#### Set-Up Git
The distributed version control system [Git](https://git-scm.com/) should be preinstalled on your computer. However, before you can use it, you should first configure some important settings.
> **Note:** If Git is not preinstalled, you can download it [here](https://git-scm.com/downloads).

###### Your Identity
The first thing you should do is to set your identity, i.e., your username and email address. Your identity is important because every Git commit you create uses this information, and once a commit is finished this information is unchangeable.
Please set your username to your full name, i.e., first name followed by last name,
and your email address to your institute email address with the following commands:

```bash
git config --global user.name "Max Mustermann"
git config --global user.email mustermann@lnm.mw.tum.de
```
> **Note:** You may want to use a different name or email address for other projects your are working on. For that purpose, you can run the above commands without the `--global` option when you are in a project folder.

[↑ Contents](#contents)

###### Default Text Editor
You can configure the default text editor that will be used whenever you need to write a message in Git.
The following command will set your default text editor to `kwrite`, a gui-based editor.

```bash
git config --global core.editor kwrite
```

> **Note:** Another popular choice is `vim`.

If you choose not to set a editor, Git will use your system’s default editor.

[↑ Contents](#contents)

###### Check Settings
To confirm the correct setup of Git, you may check your configuration settings with:

```bash
git config --list
```

[↑ Contents](#contents)

#### Set-Up GitLab
[GitLab](https://gitlab.lrz.de/) is a web-based service to manage Git repositories.
In addition to hosting the actual repositories, GitLab provides wikis, issue tracking, and an easy and transparent way for code review via [merge requests](https://gitlab.lrz.de/baci/baci/blob/master/CONTRIBUTING.md#merging-changes-into-master).
The GitLab we are using is hosted by the Leibniz-Rechenzentrum ([LRZ](https://www.lrz.de/)).
Before your start working in the BACI GitLab repository, please set up your account according to the description below.

###### User Profile
Your [profile settings](https://gitlab.lrz.de/help/user/profile/index.md) are available from the up-right corner menu bar (look for the user's avatar) > 'Settings', or under https://gitlab.lrz.de/profile.

1. Under 'Main Settings' please enter your full name, i.e., first name followed by last name, in the 'Name' field.
1. You may want to add a foto of you as a 'Public Avatar' so people can recognize you more easily.
1. Hit 'Update profile settings'.

[↑ Contents](#contents)

###### Change Your Username
Please set your `username` to something recognizable. It is recommended to set it to: first letter of first name followed by last name, all lowercase, e.g., Max Mustermann &rightarrow; `mmustermann`.

> **Note:** Your `username` is a unique namespace related to your user ID. Changing it can have unintended side effects. If you have already been using LRZ GitLab, read [how redirects will behave](https://gitlab.lrz.de/help/user/project/index.md#redirects-when-changing-repository-paths) before proceeding.

To change your `username`:

1. Navigate to your profile's 'Settings' > 'Account', or go to https://gitlab.lrz.de/profile/account.
1. Enter a new username under 'Change username'.
1. Hit 'Update username'.

[↑ Contents](#contents)

###### Notification Emails
To set your institute email address as your GitLab notification email:
1. Go to profile's 'Settings' > 'Emails' or to https://gitlab.lrz.de/profile/emails .
1. Enter your email address in the 'Email' field.
1. Hit 'Add email address'.
1. Go to profile's 'Settings' > 'Notifications', or to https://gitlab.lrz.de/profile/notifications .
1. In the 'Notification email' dop down menu choose your preferred notification email address.

To change your notification settings:
1. Go to profile's 'Settings' > 'Notifications', or to https://gitlab.lrz.de/profile/notifications .
1. Adjust your 'Global notification level' according to your preferences.  It is recommended to set your notification level at least to 'On mention'.
1. Alternatively, you may also adjust the notification level for each of your 'Groups' or 'Projects' individually.

[↑ Contents](#contents)

###### SSH Keys
[SSH keys](https://gitlab.lrz.de/help/ssh/README) allow you to establish an easy and secure connection between your computer and GitLab to push your local changes to the LRZ GitLab server.

To add a SSH key to your GitLab account please follow the instructions below (or go to the [GitLab documentation](https://gitlab.lrz.de/help/ssh/README)):

_1. Check for an existing SSH key pair_

Run the following command to check for an existing SSH key pair:

```bash
cat ~/.ssh/id_rsa.pub
```

If you see a string starting with `ssh-rsa` you already have a SSH key pair.
You should skip the next step and go directly to the copy to clipboard step.

_2. Generate a new SSH key pair_

To generate a new SSH key pair, execute the following command:

```bash
ssh-keygen -t rsa -b 4096-f ~/.ssh/id_rsa
```
You will be prompted to input a password to secure your SSH key pair. You can skip creating a password by pressing enter.
> **Note:** It is best practice to use a password, but it is not required.

_3. Copy your public SSH key to the clipboard_

Repeat step one. You should now see your public SSH key: a string starting with `ssh-rsa`.
Highlight the string and press <kbd>Ctrl</kbd> + <kbd>Shift</kbd> + <kbd>C</kbd> to copy it from the console to your clipboard.

_4. Add Your Public SSH Key to GitLab_

Navigate to profile's 'Settings' > 'SSH Keys', or go to https://gitlab.lrz.de/profile/keys.
Paste your key in the 'Key' section by pressing <kbd>Ctrl</kbd> + <kbd>v</kbd>.
Give it a relevant 'Title' and hit 'Add key'.

_5. Test your setup_

To test wether you have added your SSH key correctly, run the following command:
```bash
ssh -T git@gitlab.lrz.de
```
You should see a `Welcome to GitLab` message.

[↑ Contents](#contents)

## License

ADD A LICENSE!!

[↑ Contents](#contents)
