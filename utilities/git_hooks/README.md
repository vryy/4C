# GIT-Hooks

## Contents
1. [General background](#General background)
   1. [Set-up Git](#set-up-git)
   1. [Set-up LRZ GitLab](#set-up-lrz-gitlab)
   1. [Clone the Repository](#clone-the-repository)
   1. [Set Up the Environment](#set-up-the-environment)
   1. [Configure and Build](#configure-and-build)
   1. [Updating BACI](#updating-baci)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing](#contributing)
1. [License](#license)


## General background
From https://git-scm.com/book/pl/v2/Customizing-Git-Git-Hooks:
There are two groups of these hooks: client-side and server-side.
Client-side hooks are triggered by operations such as committing and merging, while server-side hooks run on network operations such as receiving pushed commits. 
You can use these hooks for all sorts of reasons.

### Committing-Workflow Hooks
From: https://git-scm.com/book/pl/v2/Customizing-Git-Git-Hooks:

1. The pre-commit hook: First hook that is run (even before typing commit message). Exiting non-zero from this hook aborts the commit, although you can bypass it with git commit --no-verify.Suitable for style checks e.g. trailing whitespace o even for appropriate documentation on new methods.
1. The pre-commit-msg: Edits the default message before the commit author sees it. Generally not useful for normal commits. Good for commits where the default message is auto-generated, such as templated commit messages, merge commits, squashed commits, and amended commits.
1. The commit-msg: Takes one parameter as input which is path to a temporary file that contains the commit message written by the developer. If this script exits non-zero, Git aborts the commit process.
 Useful to validate project state or commit message before allowing a commit to go through.
1. The post-commit:



## clang

## Possible modifications and variants

##
