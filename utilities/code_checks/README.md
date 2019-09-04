# GIT-Hooks

## Contents
1. [LNM previous code checks](#LNM previous code checks)
1. [New work-flow in Git: Git-hooks](#New work-flow in Git: Git-hooks)
  1. [Client-side hooks](#Client-side hooks)
  1. [Server-side hooks](#Server-side hooks)
1. [Clang-format for code-style checks](#Clang-format for code-style checks)
  1. [Useful commands](#Useful commands)
  1. [Integration in personal IDE](#Integration in personal IDE)
1. [Git-hook configuration in repository](#Git-hook configuration repository)
1. [Setup in GitLab](#Setup in GitLab)



## LNM previous code checks
Before BACI was moved to GitLab, respectively used Git for version control, a pre-commit and post-commit script was used within SVN to call several Python scripts that:
* Checked for trailing white-spaces
* Checked for tabs
* Made sure that the author of the commit had the rights to perform the commit
* Made sure that C/C++ files are dos line-endings-free
* Checked for compliance with pre-defined header format
* Checked the commit messages
* Check which files are ignored from commits and changes
* Organized the e-mails and notifications for contributors and users

As the available scripts were written for SVN specific syntax and Git enables a cleaner organization of code checks, a new routine for Git was implemented.
The latter is controlled by so called *Git-hooks* scripts which run in different stages of the Git work-flow.

## New work-flow in Git: Git-hooks
From: https://git-scm.com/book/pl/v2/Customizing-Git-Git-Hooks  <br>In general one can distinguish between two classes of hooks: <br>
**Client-side hooks** and **server-side hooks**.
In the current implementation we only use client-side hooks to control the commit behavior. This guarantees a direct feedback to the contributor and prevents too large and confusing code corrections.
Hooks can be used to run all sorts of checks before or after code is committed, merged or pushed. A short overview is given in the following:

### Client-side hooks
Client-side hooks are used on the local repository to control the commit behavior. There are different types of client-side hooks existent:
1. The **pre-commit** hook: First hook that is run (even before typing commit message). Exiting non-zero from this hook aborts the commit, although you can bypass it with `git commit --no-verify`. Suitable for style checks e.g. trailing white-space or even for appropriate documentation on new methods.
1. The **pre-commit-msg**: Edits the default message before the commit author sees it. Generally not useful for normal commits. Good for commits where the default message is auto-generated, such as template commit messages, merge commits, squashed commits, and amended commits.
1. The **commit-msg**: Takes one parameter as input which is path to a temporary file that contains the commit message written by the developer. If this script exits non-zero, Git aborts the commit process.
 Useful to validate project state or commit message before allowing a commit to go through.
1. The **post-commit**: Runs after the commit procedure. Can be used for notification of special commits or similar tasks.

### Server-side hooks
Server-side hooks run before and after pushes to the server. They can be used to set up a push policy.
Currently we do not have any server-side hooks in place.
Types of server-side hooks are:
1. **Pre-receive** hook: First script executed after push. Can be used to control references and files that would be modified due to the push
1. **Update hook**: Similar to pre-receive hook but it runs once for each branch trying to be updated.
1. **Post-receive** hook: Runs after completion of push process. Can be used for e-mail lists or within the use of a ticket tracking system. Cannot be used to stop push process.


## Clang-format for code-style checks
Clang-format is an advanced software tool that is part of a compiler front end within the LLVM release cycle. It is used to format C/C++/Java/JavaScript/Objective-C/Protobuf code.
There are several pre-existing code styles available but one can also define an own style. The tool will check existing codes for compliance with the predefined style and gives the possibility to change the code, accordingly. There are plenty of powerful options for code formatting possible and the capabilities of the tool exceed previous used checks (for trailing white-spaces or tabs) by far.
The interested reader might refer to [Format Style Options](https://clang.llvm.org/docs/ClangFormatStyleOptions.html) for more details on different style options.

### Useful commands
From https://clang.llvm.org/docs/ClangFormat.html:
<br>Clang is executed via the command `<path/to/clang/>clang-format [options] [<file> ...]`.
Useful `options` and their meaning are listed in the following:
* `-i`: files that are not compliant with the current clang format are edited in place
* `-style=<string>`: Specifies different coding styles that should be used by clang. Options for `<string>` are:
  * `file`: This option reads in a custom formatting style that is saved in a file named `.clang-format` that has to be located in the top-level folder of the project.
  * `llvm`
  * `google`: This is the style that Google uses for their projects (We use a slight modification of that style).
  * `chromium`
  * `mozilla`
  * `webkit`
* `-dump-config`: Dump the specified configuration to stdout. This command should be used with the `-style` option for the desired output style that can be saved to a `.clang-format` file. To write out a `.clang-format` file with the Google style format one can use the command `clang-format -style=google -dump-config > .clang-format`.
* `-verbose`: shows the list of processed files

It is possible to shield certain parts of a code from the formatting by clang-format. The code between `// clang-format off` and `// clang format on` will not be formatted.

### Integration in personal IDE
Taken from: http://clang.llvm.org/docs/ClangFormat.html#vim-integration
<br>Clang-format can also be used interactively during coding in the personal editor or IDE to directly enforce better code readability on the go.
Below are some configurations and settings for common editors:

* VIM: Run the clang-format standalone tool on your current buffer, optionally selecting regions to reformat. The integration has the form of a python-file which can be found under `clang/tools/clang-format/clang-format.py`. Clang is integrated by adding
```bash
map <C-K> :pyf <path-to-this-file>/clang-format.py<cr>
imap <C-K> <c-o>:pyf <path-to-this-file>/clang-format.py<cr>
```
to the *.vimrc* file. With this integration you can press the bound key and clang-format will format the current line in NORMAL and INSERT mode or the selected region in VISUAL mode. The line or region is extended to the next bigger syntactic entity. Alternatively, the configuration below enables clang-format only after saving of the file:
```bash
function! Formatonsave()
  let l:formatdiff = 1
  pyf ~/llvm/tools/clang/tools/clang-format/clang-format.py
endfunction
autocmd BufWritePre *.h,*.cc,*.cpp call Formatonsave()
```
* Atom: After installing clang-format, install necessary packages via shell command `apm install formatter` and `apm install formatter-clang-format` on Linux or use `brew` instead of `apm` on OSX. Specify the right *.clang-format* path in the program configuration.

* Eclipse: LLVM is currently working on an official plugin. Meanwhile the [GitHub solution CPPStyle by wangyz](https://github.com/wangzw/CppStyle) can be used. Also, take care of specifying the right *.clang-format* file in the configurations. Additionally, an **eclipse_style.xml** template is provided in `utilities/code_style`

* Visual Studio: Download the latest Visual Studio extension from the [alpha build site](http://llvm.org/builds/). The default key-binding is `Ctrl-R`,`Ctrl-F`.

## Git-hook configuration in repository
Currently two different client-side hooks are used in the work-flow. A *pre-commit* hook and a *commit-msg* hook. See [Client-side hooks](Client-side hooks) for reference. Both files can be found in `utilities/code_checks`. It was decided to check the code style in the pre-commit hook that consists of two parts: The clang-format style changes and a subsequent check for correct header format. After the user types `git commit` the hook gets activated. First clang-format is called and the code style formation is done automatically. The program clang format is shared in the repository and is located in `utilities/code_checks`. The associated format file *.clang-format* is located in the top level folder. For the code style the *Google-style template* is used with some slight modifications that can be seen when opening the file itself.

Afterwards the user gets a short terminal output of the changes done by clang format. If the code was already compliant to the defined clang style no extra output is given. For the header checks three different Python files namely *baciheader.py*, *header-check.py* and *inputheader.py* are called. Headers are checked for cpp-files and input-files, separately:
* cpp-files:
  * file tag
  * maintainer tag
  * level tag
* input files:
  * check .dat-files for header
  * check for maintainer

## Setup in GitLab
Taken from, [GitLab website](https://docs.gitlab.com/ee/administration/custom_hooks.html):
Standard Git hooks are located in project `hook` directory and GitLab automatically creates a symlink. Custom hooks are implemented a little differently but behave the same:
1. On the GitLab server navigate to the project's repository directory. For an installation from source the path is usually `/home/git/repositories/<group>/<project>.git`. For Omnibus installs the path is usually `/var/opt/gitlab/git-data/repositories/<group>/<project>.git`.
1. Create a new directory in this location called `custom_hooks`.
1. Inside the new `custom_hooks` directory, create a *file with a name matching the hook type*. For a pre-receive hook the file name should be pre-receive with no extension.
1. Make the hook file executable and make sure it's owned by git.
1. Write the code to make the Git hook function as expected. Hooks can be in any language. Ensure the 'shebang' at the top properly reflects the language type. For example, if the script is in Ruby the shebang will probably be `#!/usr/bin/env ruby`

Custom Git hooks must be configured on the file-system of the GitLab server. Only GitLab server administrators will be able to complete these tasks.

## Compliant file-maintainers in Baci
Our Git-hooks check for compliant file-maintainers in source and input files and will reject commits with incompliant file-maintainters. A definition of the latter is given in the following:
- A **file-maintainer** (not to be confused with **baci_maintainers**) is the person responsible for a specific source or input file in Baci
- A file-maintainer is assigned in the header declaration of a file via a:
  - `\maintainer` tag in **source-files**, followed by the name of the file-maintainer
  - `// Maintainer:` tag in **input-files**, followed by the name of the file-maintainer

  which has to be compliant with names listed in `utilities/code_checks/baci_developers.json`
- All members of the Gitlab **baci_developers** group are compliant file-maintainers
- We only allow for one compliant file-maintainer per file
- An up-to-date `baci_developers.json` list can be created by the following steps:
    1. Create a personal access token for your Baci Gitlab repository as described in the [gitlab documentation](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html)
    1. Download the list of all members of the **baci_developers** group using the command
`curl -H "Private-Token: <Token>" "https://gitlab.lrz.de/api/v4/groups/13552/members/all?per_page=100&page=1" > baci_developers.json`
    1. Commit and merge the changes in `utilities/code_checks/baci_developers.json`
