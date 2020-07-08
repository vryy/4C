#!/usr/bin/env python
import os
import sys
import argparse
import shlex
import json

from utils import get_mpi_include_path


def do_vscode_configuration(src_dir, build_dir, cmake_config):
    print("Configure VS Code")
    vscode_dir = os.path.join(src_dir, ".vscode")
    if not os.path.isdir(vscode_dir):
        # the user does not use VS Code
        return

    cmake_args = shlex.split(cmake_config)

    compiler_path = None
    compile_commands = os.path.join(build_dir, "compile_commands.json")

    # read compiler
    for cmake_arg in cmake_args:
        if cmake_arg.startswith("CMAKE_CXX_COMPILER:FILEPATH="):
            compiler_path = cmake_arg.split("=")[1]

    # check if compiler was found
    if compiler_path is None:
        print(
            "Compiler not found in the cmake configuration string. Cannot setup VS Code automatically."
        )
        return

    # check, whether the mpi header directory exists
    mpi_header = get_mpi_include_path()

    if mpi_header is None:
        # Could not find openmpi (the location is currently hard-coded)
        print(
            "The mpi include dir could not be found for your system. Please add the include dir to the possible paths in utils.py"
        )
        return

    # Check, whether there is already a configuration or not
    if not os.path.isfile(os.path.join(vscode_dir, "c_cpp_properties.json")):
        # start from scratch
        cpp_properties = {"configurations": [], "version": 4}
    else:
        # use existing configuration
        with open(os.path.join(vscode_dir, "c_cpp_properties.json"), "r") as file:
            cpp_properties = json.load(file)

        # check for compatibility
        if "version" not in cpp_properties:
            print(
                "The c_cpp_properties.json file is corrupt! Have to abort automaic setup of VS Code"
            )
            return
        elif cpp_properties["version"] != 4:
            print(
                "Currently only version 4 of c_cpp_properties.json is supported. You may need to adapt adapt_vscode.py to support more versions."
            )
            return

    # check, whether cmake wrote the compilation database
    if not os.path.isfile(compile_commands):
        print("CMAKE has not written the compilation database (compile_commands.json).")
        return

    name = os.path.basename(build_dir)

    # check whether the target was already configured
    for configuration in cpp_properties["configurations"]:
        if configuration["name"] == name:
            print("VS Code is already configured")
            return

    # create new configuration
    new_configuration = {
        "name": name,
        "includePath": ["${{workspaceFolder}}/**", mpi_header],
        "compileCommands": compile_commands,
        "defines": [],
        "compilerPath": compiler_path,
        "cStandard": "c11",
        "cppStandard": "c++11",
        "intelliSenseMode": "gcc-x64",
    }

    # add configuration
    cpp_properties["configurations"].append(new_configuration)

    # write configuration file
    with open(os.path.join(vscode_dir, "c_cpp_properties.json"), "w") as file:
        json.dump(cpp_properties, file, sort_keys=True, indent=4)

    print("Done configuring VS Code")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Setup of VS Code to be used with baci"
    )
    parser.add_argument("src_dir", help="Path to the source directory"),
    parser.add_argument("build_dir", help="Path to the build directory"),
    parser.add_argument("cmake_config", help="CMake configuration.")

    args = parser.parse_args()

    do_vscode_configuration(args.src_dir, args.build_dir, args.cmake_config)
