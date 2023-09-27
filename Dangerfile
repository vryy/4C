### General
warn("Please assign yourself (and any other authors) to this MR.", sticky: true) unless gitlab.mr_json["assignee"]

### Code owners
git_affected_files = (git.modified_files + git.added_files + git.deleted_files + git.renamed_files).join(" ")

# Do not mention @baci/baci_developers as this would create too much noise. Anyone can approve such a rule.
codeowners = `python3 ./utilities/danger/all_code_owners.py .gitlab/CODEOWNERS --files #{git_affected_files} --exclude_owners @baci/baci_developers`
message("Mentioning the code owners of files affected by this MR: #{codeowners}")

markdown("_Note: This comment will be updated when the job or pipeline is run again._")
