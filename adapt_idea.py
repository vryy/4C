#!/usr/bin/env python
import sys
import os.path

try:
    from lxml import etree
except ImportError:
    print "\n Error: python-lxml is not installed. For installation type as root:"
    print "   yum install python-lxml "
    print " exiting now...\n"
    sys.exit(1)

def adapt(cmake_cmd_line,build_type,path_to_settings):
    """update CLion's .idea/workspace.xml settings if existing"""

    if not os.path.isfile(path_to_settings):
        # user does not seem to use CLion
        return

    with open(path_to_settings,"r") as f:
        project = etree.fromstring(f.read())

    found = False
    for configurations in project.iter("configurations"):
        for currentConfiguration in configurations.iter("configuration"):
            if currentConfiguration.get("PROFILE_NAME").lower() == build_type.lower():
                print "Configuring CLion for "+build_type
                found = True
                currentConfiguration.set("GENERATION_OPTIONS", cmake_cmd_line)
    
    if not found:
        print "Warning: Could not find build config "+build_type+" in "+path_to_settings
        return

    with open(path_to_settings, "w") as fo:
        fo.write(etree.tostring(project, encoding="UTF-8"))

    print "++ Update of .idea/workspace.xml file done"

if __name__=='__main__':
    # argument order: cmake_cmd_line, buildtype, settings file
    adapt(sys.argv[1],sys.argv[2],sys.argv[3])
