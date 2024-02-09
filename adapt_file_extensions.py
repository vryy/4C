import re, os, sys
import fileinput


""" 
Re to find all baci include headers in a line
"""


def find_baci_header_includes(line):

    found = re.search('#include\s+"baci_', line)
    if found:
        return True
    else:
        return False


"""
Re to find the header guard lines define and ifndef in a line(string)
"""


def find_baci_header_guards(line):

    ifndef = re.search("#ifndef\s+BACI_", line)
    define = re.search("#define\s+BACI_", line)
    if define or ifndef:
        return True
    else:
        return False


""" 
Apply replacement on line, if the condition is statisfied
"""


def rename_in_file(inputfile, condition, pattern, replacement):

    # search for lines to replace in file
    for line in fileinput.input(inputfile, inplace=1):
        if condition(line):
            line = line.replace(pattern, replacement)
        sys.stdout.write(line)


"""
loop through folder and subfolders to update the extensions in a file
Parameters
----------
root: path
    path to starting directory
file_endling_list: list
    list containing all files which should be modified
replace_extensions: dict
    key value pair containing the old extension as key and the new extension as value

file_endling_list - list containing all files which should be modified
replace_extensions - key value pair containing the old extension as key and the new extension as value
"""


def recursive_in_file_renames(root, file_ending_list, replace_extension):

    # loop through file system
    for root, dirs, files in os.walk(root, topdown=False):

        # check if file ends with
        for name in files:
            for ending in file_ending_list:
                if name.endswith(ending):
                    for (
                        replace_extension_key,
                        replace_extension_value,
                    ) in replace_extension.items():
                        rename_in_file(
                            os.path.join(
                                root,
                                name.replace(
                                    replace_extension_key, replace_extension_value
                                ),
                            ),
                            find_baci_header_includes,
                            replace_extension_key,
                            replace_extension_value,
                        )

        for name in dirs:
            recursive_in_file_renames(
                os.path.join(root, name), file_ending_list, replace_extension
            )


"""
rename a file according to the new extension (and removes automatically the old file)
"""


def rename_file(file, new_extension):

    # get old file name
    old_name = os.path.abspath(file)

    # get file name without extension
    only_name = os.path.splitext(file)[0]

    # construct full file path and new name
    new_name = os.path.join(os.path.dirname(file), only_name + new_extension)

    # renaming the file
    os.rename(old_name, new_name)


"""
Recursive rename of all files with extension and update the Header guards accordingly.

Parameters
----------
root: path
    path-to-directory
file_extensions: dict
    key value pair containing the old extension as key and the new extension as value
header_guards_in_file: dict
    key value pair containing the old header guard as key and the new header guard as value
"""


def recursive_rename(root, file_extensions, header_guards_in_file):

    # loop through directory
    for root, dirs, files in os.walk(root, topdown=False):

        # file loop
        for name in files:
            for file_extension_key, file_extension_value in file_extensions.items():

                # check if we need to update the file
                if name.endswith(file_extension_key):

                    # create new file
                    rename_file(os.path.join(root, name), file_extension_value)

                    # update the header guards of the new file
                    for (
                        header_guard_key,
                        header_guard_value,
                    ) in header_guards_in_file.items():
                        rename_in_file(
                            os.path.join(
                                root,
                                name.replace(file_extension_key, file_extension_value),
                            ),
                            find_baci_header_guards,
                            header_guard_key,
                            header_guard_value,
                        )

        # recursive call for remaining directorys
        for name in dirs:
            recursive_rename(
                os.path.join(root, name), file_extensions, header_guards_in_file
            )


if __name__ == "__main__":

    # dictionary with old and new expression
    replace_file_extension = {".H": ".hpp"}
    header_guards_in_file = {"_H": "_HPP"}

    # file endings which must be opened and searched due to changed extensions
    file_search_extensions = [".H", ".h", ".hpp", ".cpp"]

    # check if starting directory is provided
    if len(sys.argv) < 2:
        raise TypeError(
            "Please provide a starting directory of the recursive renaming.\n"
        )

    # starting directory
    root_dir = sys.argv[1]

    # start with renaming all files according to extension
    recursive_rename(root_dir, replace_file_extension, header_guards_in_file)

    # start renaming file extensions in the files
    recursive_in_file_renames(root_dir, file_search_extensions, replace_file_extension)
