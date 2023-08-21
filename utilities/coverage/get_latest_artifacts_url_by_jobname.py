import gitlab
import argparse
import typing


def get_instance(gitlab_instance: str, private_token: str) -> gitlab.Gitlab:
    return gitlab.Gitlab(gitlab_instance, private_token=private_token)


def get_project(instance: gitlab.Gitlab, project: str):
    return instance.projects.get(project)


def get_job_id_of_latest_successful_job_with_name(
    project, branch: str, job_name: str
) -> typing.Optional[int]:
    for pipeline in project.pipelines.list(
        ref=branch, iterator=True, order_by="id", sort="desc"
    ):
        for job in pipeline.jobs.list(iterator=True):
            if job.name == job_name and job.status == "success":
                return job.id

    return None


def main():
    parser = argparse.ArgumentParser(
        prog="get_latest_artifacts_url_by_jobname.py",
        description="Searches all successful or running pipelines and returns the url of the artifact of the latest successful job with a name on a specific branch.",
    )

    parser.add_argument("gitlab_instance")
    parser.add_argument("private_token")
    parser.add_argument("project_id")
    parser.add_argument("branch")
    parser.add_argument("job_name")

    args = parser.parse_args()

    instance = get_instance(args.gitlab_instance, args.private_token)
    project = get_project(instance, args.project_id)
    job_id = get_job_id_of_latest_successful_job_with_name(
        project, args.branch, args.job_name
    )

    if job_id is None:
        raise RuntimeError("Could not find any successful job")

    print(f"{instance.api_url}/projects/{project.get_id()}/jobs/{job_id}/artifacts")


if __name__ == "__main__":
    main()
