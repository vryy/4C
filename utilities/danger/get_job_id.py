import gitlab
import argparse


def get_instance(gitlab_instance: str, private_token: str) -> gitlab.Gitlab:
    return gitlab.Gitlab(gitlab_instance, private_token=private_token)


def get_project(instance: gitlab.Gitlab, project: str):
    return instance.projects.get(project)


def get_job_id_by_name(job_iterator, job_name: str):
    for job in job_iterator:
        if job.name == job_name:
            return job.id

    return None


def main():
    parser = argparse.ArgumentParser(
        prog="get_job_id.py",
        description="Return the job id of a job with a specific name in a specific pipeline.",
    )

    parser.add_argument("gitlab_instance")
    parser.add_argument("private_token")
    parser.add_argument("project_id")
    parser.add_argument("pipeline_id")
    parser.add_argument("job_name")

    args = parser.parse_args()

    instance = get_instance(args.gitlab_instance, args.private_token)
    project = get_project(instance, args.project_id)

    job_id = get_job_id_by_name(
        project.pipelines.get(args.pipeline_id).jobs.list(iterator=True), args.job_name
    )

    if job_id is None:
        raise RuntimeError(f"Job not found with name {args.job_name}")

    print(job_id)


if __name__ == "__main__":
    main()
