#!/usr/bin/env python

import os
import sys
import subprocess

try:
    from lxml import etree
except ImportError:
    print "\n Error: python-lxml is not installed. For installation type as root:"
    print "   yum install python-lxml "
    print " exiting now...\n"
    sys.exit(1)

import random
import string


def getCompilerPaths():
    """Get compiler paths using 'g++ -v -E -P -dD'"""
    
    # to avoid filename clashes
    save = "".join([random.choice(string.letters) for x in xrange(30)])
    complicated_name_cpp = "dummy_" + save + ".cpp"
    complicated_name_output = "dummy_" + save + ".txt"
    
    # write path information into output file
    f = open(complicated_name_cpp,"w")
    f.write("\n")
    f.close()
    os.system('g++ -v -E -P -dD '+complicated_name_cpp+' &> '+complicated_name_output)
    os.system('rm '+complicated_name_cpp)

    # parse output file and store compiler paths in set
    f = open(complicated_name_output,"r")
    paths_are_comming=False
    comppath=set()
    for l in f.readlines():
        if l.find("#include <...> search starts here:") == 0:
            paths_are_comming=True
        else:
            if l.find("End of search list") == 0:
                paths_are_comming=False
            if paths_are_comming:
                comppath.add(l.strip())
    f.close()
    os.system('rm '+complicated_name_output)
    return comppath

def getPaths(build_folder):
    """Get all include paths"""
    
    cachefile = str(build_folder) + "/CMakeCache.txt"
    
    f = open(cachefile,"r")
    
    # get paths from do-configure file
    pathlist = set()
    for l in f.readlines():
	
        if (l.find("INCLUDE_INSTALL_DIR:FILEPATH=") > -1):
            pathlist.add(l.split("=")[1][0:-1])
        if (l.find("Trilinos_DIR:FILEPATH=") > -1):
            pathlist.add(l.split("=")[1][0:-1]+"/../../../include")
    # check whether a 64-bit machine is used
    version = subprocess.Popen(["uname","-r"], stdout=subprocess.PIPE).communicate()[0]
    if 'x86_64' in version:
    	pathlist.add("/usr/include/openmpi-x86_64")
    else:
    	pathlist.add("/usr/include/openmpi/1.2.4-gcc")
    
    # add compiler paths
    pathlist.update(getCompilerPaths())
    
    return pathlist

def getSymbols(fname,build_type):
    """Get all defines flags from do-configure file"""
    f = open(fname,"r")
    symbollist = set()
    for l in f.readlines():
        words = l.split()
        if len(words) > 1:
            if words[0] == "-D":
                flags = words[1].split("=")
                flag = flags[0]
                status = flags[1]
                #valid_flag = flag[0:2] == "D_"
                #if valid_flag == True and status == "ON":
                #    symbollist.add(flag.split(":")[0])
                if status == "ON":
                    symbollist.add(flag.split(":")[0])

    if build_type == "DEBUG":
    	symbollist.add("DEBUG")
    # flags that are set in CMakeLists.txt
    symbollist.add("PARALLEL")
    symbollist.add("PARMETIS")
    symbollist.add("QHULL")
    symbollist.add("BINIO")
    symbollist.add("HAVE_FFTW")
    symbollist.add("D_SHELL8")
    return symbollist

def getDefineValue(build_folder):
    """Get all define flags that specify a value from the CMakeFiles folder"""
    f = open(build_folder+"/CMakeFiles/drt_lib.dir/flags.make","r")
    definevaluelist = set()
    for l in f.readlines():
        for w in l.split():
            if w.startswith("-DBOOST_MAJOR_VERSION") or w.startswith("-DBOOST_MINOR_VERSION"):
                flag = w[2:] #loose the -D
                definevaluelist.add(flag)
    return definevaluelist

def adapt(do_configure_file,build_folder,build_type):
    """update .cproject file if existing"""

    if os.path.isfile(".cproject"): # if file exists
        f = open(".cproject","r")

        project = etree.fromstring(f.read())

        pathset = getPaths(build_folder)
        pathlist = [x for x in pathset]
        pathlist.sort()
        #print pathlist
        symbolset=getSymbols(do_configure_file,build_type)
        symbollist = [x for x in symbolset]
        symbollist.sort()
        #print symbollist
        definevalueset = getDefineValue(build_folder)
        definevaluelist = [x for x in definevalueset]
        definevaluelist.sort()
        #print definevaluelist

        # iterate over all entries named 'option'
        found_path = False
        found_symbol = False 
        for option in project.iter("option"):
            superClass = option.get("superClass")

            if superClass[-29:] == "compiler.option.include.paths":
                #print(etree.tostring(option, pretty_print=True))
                for entry in option:
                    option.remove(entry)
                for entry in pathlist:
                    option.append(etree.Element("listOptionValue", builtIn="false", value=entry ))
                #print(etree.tostring(option, pretty_print=True))
                found_path = True

            if superClass == "gnu.cpp.compiler.option.preprocessor.def" or superClass == "gnu.c.compiler.option.preprocessor.def.symbols":
                #print(etree.tostring(option, pretty_print=True))
                for entry in option:
                    option.remove(entry)
                for entry in symbollist:
                    option.append(etree.Element("listOptionValue", builtIn="false", value=entry ))
                for entry in definevaluelist:
                    option.append(etree.Element("listOptionValue", builtIn="false", value=entry ))
                #print(etree.tostring(option, pretty_print=True))
                found_symbol = True
                
        if not found_path:
        	print "Please add manually (Eclipse) any path to the project's include path section to create an initial entry in '.cproject'"
        if not found_symbol:
        	print "Please add manually (Eclipse) any symbol to the project's symbol section to create an initial entry in '.cproject'"

        #print(etree.tostring(root, pretty_print=True))
        fo = open(".cproject","w")
        fo.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
        fo.write('<?fileVersion 4.0.0?>')
        fo.write("")
        fo.write(etree.tostring(project, pretty_print=True))
        fo.close()

        print "++ Update of .cproject file done"

if __name__=='__main__':
    adapt(sys.argv[1],sys.argv[2],sys.argv[3])
