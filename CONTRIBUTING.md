## Contributing to BACI

Thank you for your willingness to contribute to BACI. The procedure to do so is given below. 
Do ensure you are familiar with all the information in our [README.md](https://gitlab.lrz.de/baci/baci/blob/master/README.md) file, as that is necessary for understanding what follows.

> **Note:**  By contributing to BACI, you implicitly agree to our [contributor license agreement](https://gitlab.lrz.de/baci/baci/blob/master/ContributorLicenseAgreement.md).

### Contents
1. [BACI Workflow](#baci-workflow)
1. [Create a GitLab Issue](#create-a-gitlab-issue)
1. [Work an Issue](#work-an-issue)
   1. [Update the Main Development Branch](#update-the-main-development-branch)
   1. [Create a Feature Branch](#create-a-feature-branch)
   1. [Make your Changes](#make-your-changes)
   1. [Update your Branch](#update-your-branch)
   1. [Test your Changes](#test-your-changes)
1. [Merging Changes into `master`](#merging-changes-into-master)
   1. [Create a Merge Request](#create-a-merge-request)
   1. [Feedback](#feedback)
   1. [Merging a Merge Request](#merging-a-merge-request)
1. [Final Clean-Up](#final-clean-up)

## BACI Workflow

BACI development is strongly based on the [GitHub Flow](https://guides.github.com/introduction/flow/index.html) which is a branch-based workflow involving only two types of branches: the `master` branch and `feature` branches.

The most important rules are: 
- Anything in the `master` branch is always deployable, i.e. considered stable.
- Development (and bugfixes) are carried out in `feature` branches.

To incorporate a `feature` branch into the `master` branch, BACI employs GitLab's *merge request (MR)* mechanism resulting in a *merge commit*.

[↑ Contents](#contents)

## Create a GitLab Issue

Navigate to BACI's [GitLab Issues page](https://gitlab.lrz.de/baci/baci/issues) and create a new issue.
The issue can be used for any number of things &mdash; reporting a bug, suggesting an enhancement, posing a question, etc.
On the new issue creation page, select an issue template from the drop down menu
to pre-populate the *Description* field with some text.
Follow the instructions in that template to give us as much information as you can
such that the issue can be understood by your fellow developers and can be tackled as soon as it is practicable.

Issues begin their life in the **Backlog** of our [Kanban board](https://gitlab.lrz.de/baci/baci/boards)
and then move through the board from left to right.
When enough information has been gathered and work is ready to begin on an issue, drag it to the **Ready** column.
If at any point in an issue's life it becomes blocked by something (either another BACI issue, or perhaps something external), 
move the issue card into the **Blocked** column to indicate that work can't proceed until something else is dealt with.

[↑ Contents](#contents)

## Work an Issue

When work commences on an issue, move the issue card to the **In Progress** column of our [Kanban board](https://gitlab.lrz.de/baci/baci/boards). 
Then the workflow to use is the following:

### Update the Main Development Branch

To keep your local `master` branch in `baci` up-to-date with the remote `master`:
```bash
cd <path/to/baci-source-code>
git checkout master
git pull
```

where `<path/to/baci-source-code>` is the location of your local BACI repository, i.e. the BACI source code.

You want to do this before starting work on a new feature branch.

[↑ Contents](#contents)

### Create a Feature Branch

Create a local branch off of `master` in `baci` on which to make your changes:
```bash
git checkout master
git checkout -b <branchName>
```

`<branchName>` can be whatever you like, though we have some recommendations:
*  Include the issue number in it in some way, for instance, `123-<restOfBranchName>`, or `<restOfBranchName>-123`.
*  Make the branch name descriptive; that is, avoid `fixSomeStuff`, `performanceTweaks`, and generic names along those lines.
*  To indicate your branch is intended solely for your own use, include your username in the branch name somehow, as in `<username>-<restOfBranchName>` or `<restOfBranchName>-<username>`.

[↑ Contents](#contents)

### Make your Changes

Do whatever work is necessary to address the issue you're tackling.

> **Note:** Break your work into logical, compilable commits.
Feel free to commit small chunks early and often in your local repository and 
then use `git rebase -i` to reorganize your commits before sharing.

[↑ Contents](#contents)

#### Commit Messages

* Make sure your commit messages reference the appropriate GitLab issue numbers using the `#<issueNumber>` syntax.
* The **first line** of the commit message should be a descriptive title, **limited to 50 characters**.
* This is then followed by a blank line, and then the rest of the commit message is a description of the changes, 
limited to 72 characters wide.

[↑ Contents](#contents)

#### Doxygen

BACI uses [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) to generate documentation from annotated source code.
Please see [this wiki page](https://gitlab.lrz.de/baci/baci/wiki/Doxygen) for our Doxygen guidelines.

[↑ Contents](#contents)

### Update your Branch

While working on your feature in your local `<branchName>` branch in `baci`, other commits will likely make it into the remote `master` branch.  
There are a variety of ways to incorporate these changes into your local feature branch. 
Our preferred possibility is
```bash
git checkout master 
git pull
git checkout <branchName>
git rebase master
```
though there are others that are equally valid.

It might happen, that conflicts arise during the `git rebase master` operation.
Resolve each conflict, then continue with `git rebase --continue`.

[↑ Contents](#contents)

### Test your Changes

To ensure your changes haven't broken anything, you'll want to run `ctest` in your BACI build directory.
A small set of test cases can be run via `ctest -L minimal`.

[↑ Contents](#contents)

## Merging Changes into `master`

To merge changes into `master`, a feature branch needs to satisfy these conditions:
* Passing code check, e.g. no trailing white spaces, proper Doxygen style, ...
* No build errors and passing all tests
* Passing code inspection by one of your fellow developers

### Create a Merge Request

When your changes are ready to be integrated into `baci`'s `master` branch,
move the issue card from **In Progress** to **Under Review** on our 
[Kanban board](https://gitlab.lrz.de/baci/baci/boards) and then:

*  Push your local feature branch up to the remote with `git push --set-upstream origin <branchName>`.
*  Navigate to BACI on GitLab and [create a new merge request](https://gitlab.lrz.de/baci/baci/merge_requests/new):
   * Be sure you choose:
      * source branch: `<branchName>`
      * target branch: `master`
   * On the new merge request creation page, select a merge request template from the dropdown menu to pre-populate the *Description* field with some text. Follow the instructions in that template to give us as much information as you can such that we can review and approve the merge request as soon as it is practicable.
* Trigger the execution of the test suite manually:
   * Go to BACI's [CI/CD](https://gitlab.lrz.de/baci/baci/pipelines) page
   * Select your branch `<branchName>` and start pipeline

[↑ Contents](#contents)

#### Work-in-Progress Merge Requests

If work on an issue is not yet complete, but you'd like to get another set of eyes on your work sooner rather than later,
you can create a "work-in-progress" merge request.
Simply begin the *Title* with `WIP:` and that will indicate to the team that this is ongoing work 
that is not necessarily meant for review yet.

If you are working on a feature addition that is fairly substantial (say greater than a month of work), 
consider creating WIP MRs.  These can be reviewed, but then you can close them without merging in the changes.
When work is complete, create a MR that includes all the changes, 
and mention all the sub-WIP-MRs that have already been reviewed in the *Description*.
This makes it easy for a reviewer to see that all the changes have already been reviewed along the way, 
rather than having to look at the entire change set at once.

[↑ Contents](#contents)

### Feedback

At this point you'll enter into a stage where you and various BACI developers will iterate back and forth until your changes are in an acceptable state and can be merged in. If you need to make changes to your merge request, make additional commits on your `<branchName>` branch and push them up to the remote.  Make sure you don't delete your remote feature branch before your merge request has been merged.

[↑ Contents](#contents)

### Merging a Merge Request

Before a MR is merged, we recommend using `git rebase -i` to squash your feature branch into the smallest number of logical commits.  Much of the time this will just mean a single commit, but you may wish to keep more than one &mdash; for instance, have the majority of your feature addition in one commit, but keep some performance tweaks separate in another commit, in case it becomes necessary to revert the performance tweaks later while keeping the rest of the feature.

Once the feature branch is ready to be merged, use the "Merge" button on the MR page on GitLab.

> **Note** Only master users can merge into the `master` branch.

If your MR *Description* has some form of "closes #\<issueNumber\>" in it somewhere, merging the MR will automatically close the associated issue, which will move the issue card from **Under Review** to **Done** on the [Kanban board](https://gitlab.lrz.de/baci/baci/boards). If not, you'll need to make this move manually and adapt each issues' labels manually.

[↑ Contents](#contents)

## Final Clean-Up

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