# Contributing to BACI

Thank you for your willingness to contribute to BACI.
The steps outlined in [Setup and Initial Configuration](#setup-and-initial-configuration) have to be performed only once, while the  [Baci Development Workflow](#the-baci-development-workflow) has to be cycled for every bit of code development in Baci.

> **Note:**  By contributing to BACI, you implicitly agree to our [contributor license agreement](https://gitlab.lrz.de/baci/baci/blob/master/ContributorLicenseAgreement.md).

BACI development is strongly based on the [GitHub Flow](https://guides.github.com/introduction/flow/index.html) which is a branch-based workflow involving only two types of branches: the `master` branch and `feature` branches.
The most important rules are:

- Anything in the `master` branch is always deployable, i.e. considered stable.
- Development (and bugfixes) are carried out in `feature` branches.

To incorporate a `feature` branch into the `master` branch, BACI employs GitLab's *merge request* mechanism resulting in a *merge commit*.

### Contents
1. [Setup and Initial Configuration](#setup-and-initial-configuration)
1. [The Baci Development Workflow](#the-baci-development-workflow)
   1. [Create a GitLab Issue](#create-a-gitlab-issue)
   1. [Work an Issue](#work-an-issue)
      1. [Create a Feature Branch](#create-a-feature-branch)
      1. [Make your Changes](#make-your-changes)
      1. [Integrate changes from `master` into your feature branch](#integrate-changes-from-master-into-your-feature-branch)
      1. [Test your Changes](#test-your-changes)
   1. [Merging Changes into `master`](#merging-changes-into-master)
      1. [Push your branch to GitLab](#push-your-branch-to-gitlab)
      1. [Create a Merge Request](#create-a-merge-request)
      1. [Feedback](#feedback)
      1. [Merging a Merge Request](#merging-a-merge-request)
   1. [Final Clean-Up](#final-clean-up)

## Setup and Initial Configuration


### Read the [README.md](https://gitlab.lrz.de/baci/baci/blob/master/README.md)

Do ensure you are familiar with all the information in our [README.md](https://gitlab.lrz.de/baci/baci/blob/master/README.md) file, as that is necessary for understanding what follows.
Double-check that your local Git repository is configured as described in [README.md](https://gitlab.lrz.de/baci/baci/blob/master/README.md#set-up-git).

### Setup your Integrated Development Environment

We recommend to use an Integrated Develoment Environment (IDE) for code development because it provides many convenient features and also eases to comply with our mandatory code style.
Set-up instructions for various IDEs can be found on the respective [Wiki page](https://gitlab.lrz.de/baci/baci/wikis/set-up-your-integrated-development-environment).

## The Baci Development Workflow

### Create a GitLab Issue

Navigate to BACI's [GitLab Issues page](https://gitlab.lrz.de/baci/baci/issues) and create a new issue.
The issue can be used for any number of things &mdash; reporting a bug, suggesting an enhancement, posing a question, etc.
On the new issue creation page, select an issue template from the drop down menu
to pre-populate the *Description* field with some text.
Follow the instructions in that template to give your fellow developers as much information as you can
such that the issue can be understood and tackled as soon as it is practicable.

Issues begin their life in the **Backlog** of our [Kanban board](https://gitlab.lrz.de/baci/baci/boards)
and then move through the board from left to right.
If at any point in an issue's life it becomes blocked by something (either another BACI issue, or perhaps something external),
move the issue card into the **Blocked** column to indicate that work can't proceed until something else is dealt with.
Add a comment to the blocked issue to explain why it is blocked and what needs to happen to resolve the **Blocked** status.

[↑ Contents](#contents)

### Work an Issue

When work commences on an issue, move the issue card to the **In Progress** column of our [Kanban board](https://gitlab.lrz.de/baci/baci/boards).
Then the workflow to use is the following:

#### Create a Feature Branch

> **Note:** It is important to keep your local `master` branch in `baci` up-to-date with the remote `master`. Hence, creation of a feature branch consists of *two* steps.

First, update the `master` branch:
```bash
cd <path/to/baci-source-code>
git checkout master
git pull
```

where `<path/to/baci-source-code>` is the location of your local BACI repository, i.e. the BACI source code.

Then, create a local branch off of `master` in `baci` on which to make your changes:
```bash
git branch <branchName>
git checkout <branchName>
```

The branch name `<branchName>` can be whatever you like, though we have some recommendations:
*  Include the issue number in it in some way, for instance, `123-<restOfBranchName>`, or `<restOfBranchName>-123`.
*  Make the branch name descriptive; that is, avoid `fixSomeStuff`, `performanceTweaks`, and generic names along those lines.
*  To indicate your branch is intended solely for your own use, include your username in the branch name somehow, as in `<username>-<restOfBranchName>` or `<restOfBranchName>-<username>`.

> **Note:** Use `git branch` to list all available branches.

[↑ Contents](#contents)

#### Make your Changes

* Do whatever work is necessary to address the issue you're tackling.
* Commit your changes locally.

> **Note:** Break your work into logical, compilable commits.
Feel free to commit small chunks early and often in your local repository and
then use `git rebase -i` to reorganize your commits before sharing.

[↑ Contents](#contents)

##### Commit Messages

* The **first line** of the commit message should be a descriptive title, **limited to 50 characters**.
* This is then followed by a blank line, and then the rest of the commit message is a description of the changes,
limited to 72 characters wide.
* Make sure your commit messages reference the appropriate GitLab issue numbers using the `#<issueNumber>` syntax.

Use the commit message to explain the context and reasons for your changes, i.e. the **What?** and **Why?**, rather than providing details on how you changed the code. Use comments in the source code or in the Doxygen documentation to explain *how* the code works.

[↑ Contents](#contents)

##### Doxygen

BACI uses [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) to generate documentation from annotated source code.
Please see [this wiki page](https://gitlab.lrz.de/baci/baci/wikis/Doxygen) for our Doxygen guidelines.

[↑ Contents](#contents)

#### Integrate changes from `master` into your feature branch

While working on your feature in your local `<branchName>` branch in `baci`, other commits will likely make it into the remote `master` branch.  There are a variety of ways to incorporate these changes into your local feature branch. Our preferred possibility is
```bash
git checkout master
git pull
git checkout <branchName>
git merge master
```
though there are others that are equally valid.

> **Note:** It might happen that conflicts arise during the `git merge master` operation. After seeing a conflict, you can do two things:
>
> * Decide not to merge. Run `git merge --abort` to abort the merge operation and to restore your version of the code (without incorporating changes from `master`).
> * Resolve the conflicts. Git will mark the conflicts in the working tree. Edit the files into shape and `git add` them to the index. Use `git commit` to seal the deal.

> **Note:** You might want to integrate changes from `master` on a regular basis to ease resolving of possible conflicts.

[↑ Contents](#contents)

#### Test your Changes

To ensure your changes haven't broken anything, run `ctest` in your BACI build directory.
A small set of test cases can be run via `ctest -L minimal`.

[↑ Contents](#contents)

### Merging Changes into `master`

To merge changes into `master`, a feature branch needs to satisfy these conditions:
* Passing code check, e.g. no trailing white spaces, proper Doxygen style, ...
* No build errors and passing all tests
* Passing code inspection by one of your fellow developers

#### Push your branch to GitLab

To publish your changes and make them available to others, you have to push them to GitLab. Before pushing your branch to GitLab, use interactive rebasing via `git rebase -i` to squash the commits on your feature branch into the smallest number of logical commits.  Much of the time this will just mean a single commit, but you may wish to keep more than one &mdash; for instance, have the majority of your feature addition in one commit, but keep some performance tweaks separate in another commit, in case it becomes necessary to revert the performance tweaks later while keeping the rest of the feature.

Push your local feature branch up to the remote with:
```bash
git push --set-upstream origin <branchName>
```

> **Important**: Use `git rebase -i` only on commits that haven't been pushed to the remote, yet.

#### Create a Merge Request

When your changes are ready to be integrated into `baci`'s `master` branch,
move the issue card from **In Progress** to **Under Review** on our
[Kanban board](https://gitlab.lrz.de/baci/baci/boards) and then:

*  Navigate to BACI on GitLab and [create a new merge request](https://gitlab.lrz.de/baci/baci/merge_requests/new):
   * Be sure you choose:
      * source branch: `<branchName>`
      * target branch: `master`
   * On the new merge request creation page, select a merge request template from the dropdown menu to pre-populate the *Description* field with some text. Follow the instructions in that template to give as much information as you can such that the merge request can be reviewed and accepted as soon as it is practicable.
   * To notify interested parties and possible reviewers about your merge request, @mention possible reviewers for this merge request as well as the maintainers of all the files you've touched in the *Interested Parties / Possible Reviewers* section of the *Description* field.
* Trigger the execution of the test suite manually:
   * Go to BACI's [CI/CD](https://gitlab.lrz.de/baci/baci/pipelines) page
   * Select your branch `<branchName>` and start the pipeline

[↑ Contents](#contents)

##### Work-in-Progress Merge Requests

If work on an issue is not yet complete, but you'd like to get another set of eyes on your work sooner rather than later,
you can create a "work-in-progress" merge request.
Simply begin the *Title* with `WIP:` and that will indicate to the team that this is ongoing work
that is not necessarily meant for review yet.

If you are working on a feature addition that is fairly substantial (say greater than a month of work),
consider creating a WIP merge request.  These can be reviewed, but then you can close them without merging in the changes.
When work is complete, create a merge request that includes all the changes,
and mention all the sub-WIP-merge-requests that have already been reviewed in the *Description*.
This makes it easy for a reviewer to see that all the changes have already been reviewed along the way,
rather than having to look at the entire change set at once.

[↑ Contents](#contents)

#### Feedback

At this point you'll enter into a stage where you and various BACI developers will iterate back and forth until your changes are in an acceptable state and can be merged in. If you need to make changes to your merge request, make additional commits on your `<branchName>` branch and push them up to the remote.  Make sure you don't delete your remote feature branch before your merge request has been merged.

[↑ Contents](#contents)

#### Merging a Merge Request

Once the feature branch is ready to be merged, use the "Merge" button on the merge request page on GitLab.

> **Note** Only Maintainers can merge into the `master` branch.

If your merge request *Description* has some form of "closes #\<issueNumber\>" in it somewhere, merging the merge request will automatically close the associated issue, which will move the issue card from **Under Review** to **Done** on the [Kanban board](https://gitlab.lrz.de/baci/baci/boards). If not, you'll need to make this move manually and adapt each issues' labels manually.

[↑ Contents](#contents)

### Final Clean-Up

To keep the repository clean, delete your feature branch *after* merging it into `master`.  When you merge a merge request, GitLab will give you the option to click a button to remove the source branch.  If you click this, then
```bash
git fetch --prune
```

will remove the remote tracking information from your local repository.  Alternatively, you could skip the GitLab button and use
```bash
git push origin --delete <branchName>
```

Either way is completely fine.  After that you can remove your local branch with
```bash
git branch -D <branchName>
```

[↑ Contents](#contents)
