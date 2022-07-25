#! ./utilities/baci-python-venv/bin/python
import sys
import os.path
from lxml import etree


def adapt(cmake_cmd_line, build_type, path_to_settings):
    """update CLion's .idea/workspace.xml settings if existing"""

    if not os.path.isfile(path_to_settings):
        # user does not seem to use CLion
        return

    with open(path_to_settings, "rb") as f:
        project = etree.parse(f)

    found = False
    for configurations in project.iter("configurations"):
        for currentConfiguration in configurations.iter("configuration"):
            if currentConfiguration.get("PROFILE_NAME").lower() == build_type.lower():
                print("Configuring CLion for " + build_type)
                found = True
                currentConfiguration.set("GENERATION_OPTIONS", cmake_cmd_line)

    if not found:
        print(
            "Warning: Could not find build config "
            + build_type
            + " in "
            + path_to_settings
            + "."
            "Create it manually in CLion's project settings and rerun do-configure to adapt the CMake options."
        )
        return

    with open(path_to_settings, "wb") as fo:
        fo.write(etree.tostring(project, encoding="UTF-8"))

    print("++ Update of .idea/workspace.xml file done")


if __name__ == "__main__":
    # argument order: cmake_cmd_line, buildtype, settings file
    adapt(sys.argv[1], sys.argv[2], sys.argv[3])
