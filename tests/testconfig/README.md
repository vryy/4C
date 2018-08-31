# BACI-Testing

BACI is tested via GitLab CI/CD pipelines. Tests can be started by schedulers or whenever something changes in the repository. We test different configurations at LNM and IMCS, for a pipeline to pass it has to succeed at both LNM and IMCS. This document gives a quick overview over the important steps.

## Terminology

- `job`: A series of shell commands that will be executed on a testing machine.
- `pipeline`: A collection of jobs that will be tested. The green check mark next to a commit means, that the pipeline passed. Some pipelines only contain one job.
- `gitlab-runner`: A daemon on a testing machine that collects jobs from GitLab and runs them. Each runner has tags that are used to filter out the jobs sent to that runner. A detailed documentation can be found [here](https://docs.gitlab.com/runner/configuration/advanced-configuration.html).
- `.gitlab-ci.yml`: The configuration file for GitLab testing. In BACI it is under the `tests/testconfig/.gitlab-ci.yml`. A detailed documentation can be found [here](https://docs.gitlab.com/ee/ci/yaml/).

## How and when are pipelines started?

Every time something changes in the repository (commit, merge, etc.), GitLab checks for a `.gitlab-ci.yml` file in the repository (at the new state). This file will start ONE pipeline and the jobs in that pipeline (if no jobs are to be run, no pipeline is created). In addition to that pipelines can be triggered over the GitLab web interface and by GitLab schedulers.

## Testing in BACI

### LNM and IMCS

Since the configurations at LNM and IMCS vary, there are some environment variables defined in `.gitlab-ci.yml` that are passed to the gitlab-runner (for example at IMCS the `do-configure` script needs to be executed by `bash`, this is set in the `CTEST_CONFIGURE_PREFIX` variable).

### Daily tests

Each night a full release and a full debug test are started by a GitLab Scheduler at both LNM and IMCS. If the release test passes, `doxygen` is build and copied to a folder with constant path (if the documentation is automatically loaded to a server it is better if it is always in the same location, the path is given with a variable and can be different at LNM and IMCS).  

### Minimal tests

Every time something changes in the `master` branch, a pipeline is created that performs the code checks and a minimal test. Only if the code checks are successful, BACI is build and the minimal tests are run. The code check can be performed by any gitlab-runner that picks it up, the build and minimal tests are performed at both LNM and IMCS.

### User triggered tests

Since we allow merges to `master` only for tested commits, the user has to start a pipeline on the commit when submitting a merge request. This can be done under the GitLab web interface: Goto `CI/CD - Pipelines - Run Pipeline`, then select the branch you want to test (the latest commit on that branch will be tested) and push the button `Create pipeline`. This will first perform the code checks and if they pass a full test at both LNM and IMCS.