#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Maintainer: A. Nagler
#
# Script for searching parameters which arem't used any longer

from read_ccarat_NIGHTLYTESTCASES import read_ccarat, write_ccarat
from elements    import bcdictionary, surfaces
from progress    import progress
from sets        import Set

import sys, subprocess, re

# Dictionary of unused parameters we want to keep in the code
UNUSED_MAT_TO_KEEP = []


if __name__=='__main__':
  
    if len(sys.argv) < 2:
      print "usage: %s source-path" % sys.argv[0]
      sys.exit(1)    
    
    # collect source files of src
    files_to_search = []
    
    name = sys.argv[1]
    if name[-1] == '/':
	global_src_path = name + 'src/'
    else:	
	global_src_path = name + '/' + 'src/'
    
    # search for *.H and *.cpp files in the given src path	
    baci_heads = subprocess.check_output('find ' + global_src_path + ' -name *.H', shell=True)
    baci_heads = baci_heads.split()	
    
    baci_cpp   = subprocess.check_output('find ' + global_src_path + ' -name *.cpp', shell=True)
    baci_cpp   = baci_cpp.split()
    
    # exclude drt_validmaterials.H/.cpp from search
    baci_heads.remove(global_src_path + 'drt_inpar/drt_validmaterials.H')
    baci_cpp.remove(global_src_path + 'drt_inpar/drt_validmaterials.cpp')
    
    print 'Start to search inpar materials'
    
    # Initialization of Set for unused inpar parameters. Set used to avoid double entries
    fail = Set()
		
    inpar_materials = subprocess.check_output('grep INPAR ' + global_src_path +  'drt_inpar/drt_validmaterials.cpp', shell=True )
    inpar_materials = inpar_materials.split()
    
    # Grep for value in code. If the value doesn't appear than the input parameter might be unused
    # Check first part of files and if this failes check second part
    # If both fail than the argument doesn't exists  
	    
    for inpa_mat in progress('Searching inpar materials', inpar_materials):
	
	# special treatment of lines which include tuple<int> before the parameter
	# i.e. tuple<int>(INPAR::CAVIATION::TwoWayFull
	# method is to substitute the brackets ( with an ,
	inpa_mat = re.sub(r'\(',',', inpa_mat)
	inpa_mat = (inpa_mat.split(','))
	
	# search for desired inpar value
	for im in inpa_mat:
	    im = im.strip(' \n,;()')
	    if im[:5]== 'INPAR':
		inpa_mat = im
		break	
	# else path necessary, since the INPAR parameter values could be in comments
	else:
	    continue		
			  
	try:
	    if inpa_mat in UNUSED_MAT_TO_KEEP:
		continue
	except KeyError:
	    pass
	  
	try:
	    test = subprocess.check_output( '/bin/grep ' +  inpa_mat + " " + " ".join(baci_heads), shell=True)
	except subprocess.CalledProcessError: 
	    try:
		test = subprocess.check_output( '/bin/grep ' +  inpa_mat + " " + " ".join(baci_cpp), shell=True)
	    except subprocess.CalledProcessError: 
		fail.update([inpa_mat])	    

    if not fail:
	print "Found no unused input material in code"		
    else:
	# Retype to list in order to be able to sort
	final_fail_printout = [f for f in fail]
	final_fail_printout.sort()
	print "The following input material only exist in drt_validmaterials.cpp"
	print "\n".join(final_fail_printout)
	sys.exit(1)	