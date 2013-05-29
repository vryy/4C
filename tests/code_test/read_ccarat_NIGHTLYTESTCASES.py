#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

def read_ccarat(filename):
    section_names = []
    sections = {}
    sections["header"]=[]
    
    name_re = re.compile('^--+([^-].+)$')
    name = None

    for l in file(filename):
        line = l.strip()
        match = name_re.match(line)
        #if line[:2] == '//':
            #pass
        if match:
            name = match.group(1)
            section_names.append(name)
            #internal = name.lower().replace(' ', '_')
            sections[name] = []
        elif name:
            sections[name].append(line.split())
        ## added this to include all the other comments in the nightly test cases
        else:
            sections["header"].append(line.split())

    return section_names, sections

def write_ccarat(filename, section_names, sections):
    f = file(filename, "w")
    for name in ["header"]:
        sec = sections[name]
        for line in sec:
            #print line
            print >>f, " ".join(line)
    for name in section_names:
        if  name == "header":
            continue
        else:
            print >>f, "-"*(76-len(name)) + name
            sec = sections[name]
            for line in sec:
                if len(line)==2 and line[0] not in ["END","LAYER:","Non"]:
                    print >>f,line[0]," "*(30-len(line[0])),line[1]
                elif len(line)>3 and line[2]=='//':
                    print >>f,line[0]," "*(30-len(line[0])),line[1]," "*(30-len(line[1]))," ".join(line[2:])
                else:
                    print >>f, " ".join(line)
    #print >>f, "// END" # because it doesn't ignore the comments anymore
    f.close()

def getparam(section,name):
    for line in section:
        if line[0]==name:
            return line[1]
    return None

def setparam(section,name,value):
    found = False
    for line in section:
        if line[0]==name:
            line[1] = value
            found = True
    if not found:
        section.append([name,value])
            
if __name__ == '__main__':
    import sys
    import pprint

    if len(sys.argv) != 2:
        print "usage: %s ccarat-file" % sys.argv[0]
        sys.exit(1)

    section_names, sections = read_ccarat(sys.argv[1])
    
    #pprint.pprint(sections)
    #pprint.pprint(section_names)

    write_ccarat(sys.argv[1] + ".cpy", section_names, sections)
