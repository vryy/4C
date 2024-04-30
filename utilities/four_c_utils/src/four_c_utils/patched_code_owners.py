"Patched version of the codeowners module"

import re
from codeowners import CodeOwners


def patched_code_owners(file):
    owner = CodeOwners(file)

    for i, (pattern, path, owners, line_number, section_name) in enumerate(owner.paths):
        if pattern.pattern.startswith(r"\A"):
            if path.startswith("/"):
                # path with anchor is given in codeowners file -> only match files at exact location
                new_pattern = re.compile(pattern.pattern.replace(r"\A", r"(?:\A|\A/)"))
            else:
                # path without anchor is given in codeowners file -> also match files in subfolder
                new_pattern = re.compile(pattern.pattern.replace(r"\A", r"(?:\A|/)"))

            new_entry = (new_pattern, path, owners, line_number, section_name)
            owner.paths[i] = new_entry  # replace old entry with fixed version

    return owner
