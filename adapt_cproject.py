#!/usr/bin/env python

from lxml import etree
from sets import Set
import os
import sys

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
    comppath=Set()
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

def getPaths(fname):
    """Get all include paths"""
    f = open(fname,"r")
    
    # get paths from do-configure file
    pathlist = Set()
    for l in f.readlines():
        found = l.find("INCLUDE_INSTALL_DIR=") == 0 \
             or l.find("Trilinos_DIR=") == 0
        if (l.find("INCLUDE_INSTALL_DIR=") == 0):
            pathlist.add(l.split("=")[1][1:-2])
        if (l.find("Trilinos_DIR=") == 0):
            pathlist.add(l.split("=")[1][1:-2]+"/include")
    
    # add compiler paths
    pathlist.union_update(getCompilerPaths())
    
    return pathlist

def getSymbols(fname):
    """Get all defines flags from do-configure file"""
    f = open(fname,"r")
    symbollist = Set()
    for l in f.readlines():
        words = l.split()
        if len(words) > 1:
            if words[0] == "-D":
                flags = words[1].split("=")
                flag = flags[0]
                status = flags[1]
                valid_flag = flag[0:2] == "D_" \
                          or flag[0:5] == "DEBUG"
                if valid_flag == True and status == "ON":
                    symbollist.add(flag.split(":")[0])
    symbollist.add("CCADISCRET")
    return symbollist

def adapt(do_configure_file):
    """update .cproject file if existing"""

    if os.path.isfile(".cproject"): # if file exists
        f = open(".cproject","r")

        project = etree.fromstring(f.read())

        pathset = getPaths(do_configure_file)
        pathlist = [x for x in pathset]
        pathlist.sort()
        #print pathlist
        symbolset=getSymbols(do_configure_file)
        symbollist = [x for x in symbolset]
        symbollist.sort()

        # iterate over all entries named 'option'
        for option in project.iter("option"):
            superClass = option.get("superClass")

            if superClass[-29:] == "compiler.option.include.paths":
                #print(etree.tostring(option, pretty_print=True))
                for entry in option:
                    option.remove(entry)
                for entry in pathlist:
                    option.append(etree.Element("listOptionValue", builtIn="false", value=entry ))
                #print(etree.tostring(option, pretty_print=True))

            if superClass == "gnu.cpp.compiler.option.preprocessor.def" or superClass == "gnu.c.compiler.option.preprocessor.def.symbols":
                #print(etree.tostring(option, pretty_print=True))
                for entry in option:
                    option.remove(entry)
                for entry in symbollist:
                    option.append(etree.Element("listOptionValue", builtIn="false", value=entry ))
                #print(etree.tostring(option, pretty_print=True))

        #print(etree.tostring(root, pretty_print=True))
        fo = open(".cproject","w")
        fo.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
        fo.write('<?fileVersion 4.0.0?>')
        fo.write("")
        fo.write(etree.tostring(project, pretty_print=True))
        fo.close()

        print "++ Update of .cproject file done"

if __name__=='__main__':
    adapt(sys.argv[1])
