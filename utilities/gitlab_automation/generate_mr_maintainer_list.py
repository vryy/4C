"""
Control-module for the automatic generation of a list of maintainers for BACI GitLab Merge Requests.

The script:
1. collects the data of all open Merge Requests (MR) of the GitLab BACI Repo
2. identifies a list of the maintainers of all changed files for each open MR
3. updates the maintainer list in the description of the Merge Request accordingly
"""

import os
import json
import re
import requests
import subprocess
import sys

####################################################################################################
####################################################################################################
# Setup the script:
# 1. Add the ssh public key of the local user (on the local machine) to BACI-Testbot account.
#    This is needed to enabling passwordless cloning of the repository.
# 2. Create an access token for BACI-Testbot and save it to an empty text file on the local machine.
#    This access token will be needed to authorize the GET and PUT requests to the GitLab url.
#    It also means that all changes in the Merge Request descriptions will be attributed to the
#    BACI-Testbot account.
# 3. On the local machine, create a folder with a local BACI clone (including this file).
# 4. Setup a suitable Python environment (conda or virtualenv) on local machine.
#    Script has been tested with:
#    Python 3.7.5
#    Requests 2.22.0
# 4. Setup a cronjob on the local machine to schedule the execution of this script.
#    In the cronjob you need to supply the path to the access token as input to this script.
####################################################################################################
####################################################################################################


####################################################################################################
####################################################################################################
# !!!!!!!!!!!!!!!!!!!!! No need to touch the code below !!!!!!!!!!!!!!!!!!!!!
####################################################################################################
####################################################################################################
# based on this script's path get the base BACI directory
my_dir = os.path.dirname(os.path.abspath(__file__))
BACI_DIR = os.path.abspath(os.path.join(my_dir, '..', '..'))


MaintainerRegex = re.compile(r"[mM]aintainer:?\s*(.*)")

MaintainerList_Identifier_Top = r"<!-- begin maintainer list -->"
MaintainerList_Identifier_Bot = r"<!-- end maintainer list -->"
MaintainerListTopRegex = re.compile(MaintainerList_Identifier_Top, re.MULTILINE)
MaintainerListBotRegex = re.compile(MaintainerList_Identifier_Bot, re.MULTILINE)

# GitLab Project ID (of BACI)
BACI_GITLAB_PROJECT_ID = 22865
# GitLab url to clone BACI repository with ssh
REPO_SSH_URL = "git@gitlab.lrz.de:baci/baci.git"


def run_subprocess(cmd_list, cwd=None):
    """ Use subprocess to run a command. """
    cmd = " ".join(cmd_list)
    print("Executing:")
    print(f"    {cmd}")
    process = subprocess.Popen(
        cmd_list, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
    )
    stdout, stderr = process.communicate()
    return stdout, stderr, cmd


def get_open_mrs(access_token):
    """ Get data from all open merge request of BACI GitLab repo. """

    url = "https://gitlab.lrz.de/api/v4/projects/" + str(BACI_GITLAB_PROJECT_ID) + "/merge_requests"
    header = {"PRIVATE-TOKEN": access_token}
    data = {"state": "opened"}

    response = requests.get(url, data=data, headers=header)

    mr_data = json.loads(response.text)

    return mr_data


def get_mr_data(mr_iid, access_token):
    """ Get the Merge Request data from a specific Merge Request of BACI GitLab repo. """

    url = (
        "https://gitlab.lrz.de/api/v4/projects/"
        + str(BACI_GITLAB_PROJECT_ID)
        + "/merge_requests/"
        + str(mr_iid)
    )
    header = {"PRIVATE-TOKEN": access_token}

    response = requests.get(url, headers=header)

    mr_data = json.loads(response.text)
    error_message = mr_data.get("message", None)

    if error_message:
        raise Exception(
            f"An error occurred while contacting GitLab:\n " f"Error message:\n " f"{error_message}"
        )
    return mr_data


def prepare_repo():
    """
    Prepare the local repository for current run: includes clean up and getting latest changes.
    """

    cmd_list = ["git", "checkout", "master"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    # delete all branches except for the master branch
    cmd_list = ["git", "branch"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    branches = [line.strip(" ") for line in stdout.split("\n") if "master" not in line and line]
    for branch in branches:
        cmd_list = ["git", "branch", "-D", branch]
        stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    # clean up refs of master
    cmd_list = ["git", "fetch", "--prune", "origin"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    # get the latest changes
    cmd_list = ["git", "fetch", "--all", "origin"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    # update the master branch
    cmd_list = ["git", "reset", "--hard", "origin/master"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)


def merge_master_into_temp_work_branch(master, working, temp_work_branch):
    """ Try an automatic merge of latest master changes onto current working state. """

    # checkout a new branch where we try a to merge the latest master changes
    cmd_list = ["git", "checkout", "-b", temp_work_branch, working]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)

    if "error" in stderr or "fatal" in stderr:
        raise Exception(f"Error in:\n " f"    {cmd}\n" f"Error message:\n " f"    {stderr}")

    # try to merge into the new branch
    cmd_list = ["git", "merge", master, temp_work_branch]
    stdout, stderr, cmd = run_subprocess(cmd_list=cmd_list, cwd=BACI_DIR)

    merge_success = True
    if 'Automatic merge failed' in stdout:
        # abort the merge, if automatic merging fails
        merge_success = False
        cmd_list = ["git", "merge", "--abort"]
        stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    return temp_work_branch, merge_success


def clean_up_repo(temp_work_branch):
    """ Clean up the repository: delete temporary working branch. """

    cmd_list = ["git", "checkout", "master"]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)
    cmd_list = ["git", "branch", "-D", temp_work_branch]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)


def get_maintainer(target, temp_working_branch, merge_success=False):
    """
    Uses git diff to generate a list of changed files and gets all maintainers of the files
    """

    if BACI_DIR is None:
        raise Exception("You need to define path to your baci repo.")

    # execute git diff on both commits
    cmd_list = [
        "git",
        "diff",
        target + '...' + temp_working_branch,
        "--name-status",
        r"--find-renames=50%",
    ]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)

    print("results in:")
    print(stdout)

    if stderr:
        raise Exception(f"An error occurred during git diff.\n " f"Error message:\n {stderr}")

    # every line corresponds to a changed file, skip empty lines
    changes = [line.split("\t") for line in stdout.split("\n") if line != ""]

    extended_changes = list()
    for change in changes:
        # split all renames into delete and add
        if "R" in change[0]:
            extended_changes.append(["D", change[1], change[0]])
            extended_changes.append(["A", change[2], change[0]])
        else:
            extended_changes.append(change)

    gitlabuser_to_file = {}
    other_files = list()

    # read file maintainers
    with open(os.path.join(BACI_DIR, "utilities", "code_checks", "baci_developers.json"), "r") as f:
        maintainers_obj = json.load(f)
        maintainer_to_obj = {item["name"]: item for item in maintainers_obj}

    for change in extended_changes:
        # get maintainer of the file
        change_mode = change[0]
        original_file = change[1]

        # if file was added it can only be found in the new commit
        if "A" in change_mode:
            sources = [temp_working_branch]
        # if file was deleted it can only be found on master
        elif "D" in change_mode:
            sources = [target]
        # catch unknown change type. It should not happen! See the docs of git diff
        elif "X" in change_mode:
            raise Exception(
                'X: "unknown" change type (most probably a bug, please report it) (from git docs)'
            )
        elif "U" in change_mode:
            raise Exception(
                'Unmerged files. This script is not intended for unresolved merge commits.'
            )
        elif "T" in change_mode:
            raise Exception('Type change detected, i.e. regular file, symlink, submodule, ....')
        elif "B" in change_mode:
            raise Exception('Files with broken pairing.')
        # in all other cases check files on both master and new commit
        # M:  modification of the contents or mode of a file
        else:
            sources = [target, temp_working_branch]

        for source in sources:
            file_and_mode = [original_file, change_mode]

            maintainer = find_maintainer(original_file, source)

            if maintainer not in maintainer_to_obj:
                print(
                    "Could not find maintainer of file {0} which is {1} in the official "
                    "list".format(original_file, maintainer)
                )
                if file_and_mode not in other_files:
                    other_files.append(file_and_mode)
            else:
                maintainer_obj = maintainer_to_obj[maintainer]

                if maintainer_obj["username"] not in gitlabuser_to_file:
                    gitlabuser_to_file[maintainer_obj["username"]] = list()

                if file_and_mode not in gitlabuser_to_file[maintainer_obj["username"]]:
                    gitlabuser_to_file[maintainer_obj["username"]].append(file_and_mode)

    if merge_success:
        md_changes = "List includes latest master changes:\r\n"
    else:
        md_changes = (
            "List does NOT include latest master changes (expect conflicts during merge):\r\n"
        )

    for username, files_and_modes in gitlabuser_to_file.items():
        md_changes += "<details>\n"
        md_changes += "<summary>" + "@" + username + ":" + "</summary>\n"

        md_changes += "\n"
        for file, mode in files_and_modes:
            md_changes += "* " + " `" + file + "` (" + mode + ")\n"

        md_changes += "</details>\n"
        md_changes += "\n"

    if len(other_files) > 0:
        md_changes += "<details>\n"
        md_changes += "<summary>" + "unknown:" + "</summary>\n"

        md_changes += "\n"
        for file, mode in other_files:
            md_changes += "* " + " `" + file + "` (" + mode + ")\n"

        md_changes += "</details>\n"
        md_changes += "\n"

    return md_changes


def find_maintainer(file, source):
    """ Find the maintainer of a file in a source (i.e. a commit or branch). """
    cmd_list = ["git", "cat-file", "-p", source + ":" + file]
    stdout, stderr, cmd = run_subprocess(cmd_list, cwd=BACI_DIR)

    if stderr:
        print(stdout)
        print(stderr)
        raise Exception("An error occurred during git cat-file.")

    contents = stdout.split("\n")

    maintainer = None
    for line in contents:
        match = MaintainerRegex.search(line)
        if match is not None:
            maintainer = match.group(1).strip()
            break

    return maintainer


def update_mr_description(maintainer_list, description):
    """ Update the Merge Request description with maintainer list. """

    # get the current maintainer list
    match_maintainer_list_top = MaintainerListTopRegex.search(description)
    match_maintainer_list_bot = MaintainerListBotRegex.search(description)

    # if we can't find the identifiers, we assume there is no maintainer list
    current_maintainer_list = ""
    # current maintainer list between the identifiers
    if match_maintainer_list_top and match_maintainer_list_bot:
        idx_begin = match_maintainer_list_top.regs[0][0]
        idx_end = match_maintainer_list_bot.regs[0][1]
        current_maintainer_list = description[idx_begin:idx_end]

    # prepare the final version of the new maintainer list
    line_ending = "\r\n"
    new_maintainer_list = (
        MaintainerList_Identifier_Top
        + line_ending
        + maintainer_list
        + MaintainerList_Identifier_Bot
    )

    # change mr description only if the maintainer list has changed and there already is a list
    if current_maintainer_list != new_maintainer_list and current_maintainer_list:
        updated_description = description.replace(current_maintainer_list, new_maintainer_list)
    # there was no list yet so append the new maintainer list
    elif not current_maintainer_list:
        updated_description = (
            description + "\n\n## Maintainers of modified files:\n" + new_maintainer_list
        )
    # current and new list are identical do not do anything
    else:
        updated_description = None

    return updated_description


def put_mr_data(mr_iid, description, read_time_stamp, access_token):
    """ Update Merge Request description in GitLab repo with PUT request. """

    url = (
        "https://gitlab.lrz.de/api/v4/projects/"
        + str(BACI_GITLAB_PROJECT_ID)
        + "/merge_requests/"
        + str(mr_iid)
    )
    header = {"PRIVATE-TOKEN": access_token}
    data = {"description": description}

    # double check that the MR has not been updated since the MR data has been read
    # get the current time stamp of the last update
    current_mr_data = get_mr_data(mr_iid, access_token=access_token)
    write_time_stamp = current_mr_data["updated_at"]

    if write_time_stamp == read_time_stamp:
        response = requests.put(url, data=data, headers=header)
    else:
        raise Warning(
            "The Merge Request has been updated during the generation of list.\n"
            "Not adding the new list."
        )
        response = None

    return response


if __name__ == "__main__":

    # read the access token from externally supplied file
    if len(sys.argv) < 2:
        raise Exception("Missing input. Please supply path to access token as input.\n")
    if len(sys.argv) > 2:
        raise Warning(f"Too many input variables.\n Ignoring: {sys.arv[2:]}.\n")
    access_token_file = sys.argv[1]

    with open(access_token_file, "r") as f:
        my_access_token = f.read().rstrip()

    # get iids of all open Merge Requests
    open_mrs = get_open_mrs(access_token=my_access_token)

    for open_mr in open_mrs:
        # Update or get local clone of BACI repository
        prepare_repo()

        mr_iid = open_mr["iid"]

        # get the detailed MR data
        mr = get_mr_data(mr_iid=mr_iid, access_token=my_access_token)
        # commit to be merged into master
        head = mr["diff_refs"]["head_sha"]
        # latest commit on master
        master = mr["diff_refs"]["start_sha"]
        # shared base of master and working branch
        base = mr["diff_refs"]["base_sha"]
        # current description of the Merge Request
        mr_description = mr["description"]
        # timestamp of last update of MR
        updated_at = mr["updated_at"]

        # generate list of maintainers of modified files based on latest master changes
        temp_work_branch = "get-maintainers-mr-iid-" + str(mr_iid)
        try:
            merge_success = merge_master_into_temp_work_branch(master, head, temp_work_branch)
            # call the core routine which extracts the maintainer list
            maintainers = get_maintainer(
                target=master, temp_working_branch=temp_work_branch, merge_success=merge_success
            )
        finally:
            # delete the temporary work branch
            clean_up_repo(temp_work_branch)

        print("\n")
        print(maintainers)
        # update the maintainer list section of the Merge Request
        updated_mr_description = update_mr_description(maintainers, mr_description)
        # if the Merge Request description has been updated PUT the new description on GitLab
        if updated_mr_description:
            put_mr_data(
                mr_iid,
                updated_mr_description,
                read_time_stamp=updated_at,
                access_token=my_access_token,
            )
