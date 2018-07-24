## Contributing to BACI

This page documents at a very high level how to contribute to BACI.
We use [Git][] to develop BACI.
Git is an extremely powerful version control tool that supports many different
workflows for individual development and collaboration.
Here, we document procedures used for the development of BACI.

[Git]: http://git-scm.com

### Contents
1. [BACI Workflow](#baci-workflow)

#### BACI Workflow

BACI development is strongly based on the [GitHub Flow][] which is a branch-based workflow involving only two types of branches: the `master` branch and `feature` branches.<br/>
There's only one rule: anything in the `master` branch is always deployable, i.e. considered stable.<br/>
Preferably, a `feature` branch is used for the development of a single feature (or sometimes bug fix).

1. Anything in the `master` branch is deployable, i.e. considered stable
1. To work on something new, create a descriptively named branch off of the `master` branch with
    ```bash
     cd <BaciDir>
     git checkout master
     git pull
     git checkout -b feature-introduce-fsi
     ```
	(another feature may be developed in parallel by branching off the `master` branch again)
1. Commit to the `feature` branch locally and regularly push your work to the same named branch in the GitLab repository
1. When you think the branch is ready for merging, open a [merge request][] for your `feature` branch. You can submit a merge request through the GitLab web interface. Your request will be reviewed by someone else. If any changes are needed (due to comments in the merge request or failing tests), they are committed to the `feature` branch.
1. When the branch is ready, the feature is signed off and merged to the `master` branch.
After merging, the `feature` branch is deleted locally and on the repository.

Some special cases:
- If another merge request is accepted before you finish development on your own `feature` branch and the new commit on the `master` branch results in conflicts with your new `feature`, you have to incorporate the new commit into your `feature` branch. This can be done by [rebasing][] your `feature` branch to the new commit on the `master` branch with
```bash
 git rebase
 git push --force
 ```

 [//]: # (or by merging the `master` branch into your `feature` branch.)
- If two merge requests occur at the same time.
The default handling is: first come, first serve.
However, you may of course try to find an agreement with the other developer(s) if you find there might be a better solution for a particular situation.

Additional references:<br/>
[GitHub Flow Blog:](http://scottchacon.com/2011/08/31/github-flow.html)  first formal description in a blog post by Scott Chacon in 2011<br/>
[Git In Practice:](https://www.manning.com/books/git-in-practice) book by M. McQuaid (Chp. 14.1)

[↑ Contents](#contents)

[GitHub Flow]: https://guides.github.com/introduction/flow/index.html
[merge request]: https://gitlab.lrz.de/help/user/project/merge_requests/index.md
[rebasing]: https://git-scm.com/docs/git-rebase