#!/usr/bin/env python

import sys

def getDefineValue(build_folder):
    """Get all define flags that specify a value from the CMakeFiles folder"""
    definevalueset = set()
    
    with open(build_folder+"/CMakeFiles/drt_lib.dir/flags.make","r") as f:
      for l in f.readlines():
        if l.startswith("CXX_DEFINES"):  
          for w in l.split():
            if w.startswith("-D"):
              definevalueset.add(w[2:].replace("\\\"", "\""))
    
    return list(definevalueset)
    
def adapt(build_folder):
    """update Doxyfile.defs file"""

    definevalues = getDefineValue(build_folder)

    line = "PREDEFINED             = "
    for d in definevalues:
      line += str(d)+" "
    line += "\n"

    with open("doc/Doxyfile.defs", "w") as ff:
      ff.write(line)

    print "++ created Doxyfile definitions"

if __name__=='__main__':
    adapt(sys.argv[1])
