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

# Dictionary of unused inpar materials we want to keep in the code
UNUSED_MAT_TO_KEEP = []


if __name__=='__main__':
  
    if len(sys.argv) < 2:
      print "usage: %s source-path" % sys.argv[0]
      sys.exit(1)    
    
    # global src path
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
        
    subprocess.call('svn praise ' + global_src_path +  'drt_inpar/drt_validmaterials.cpp > praise.txt' , shell=True )	    
	    
    for inpa_mat in progress('Searching inpar materials', inpar_materials):
	
	# Exclude commentary parts
	if not inpa_mat.strip(' ')[:2] == '//':
    
	    # Get rid of leading and following brackets and of semicolons
	    inpa_mat = re.sub(r'\(',',', inpa_mat)
	    inpa_mat = re.sub(r'\)',',', inpa_mat)
	    inpa_mat = re.sub(r'\;',',', inpa_mat)
	    
	    # Search for keyword INPAR, if not found continue with loop
	    result_search = re.search('INPAR', inpa_mat)
	    if result_search:
		inpa_mat = inpa_mat[result_search.start():]
		inpa_mat = inpa_mat.split(',')[0]
	    else:
	        continue
	
	# else path necessary, since the INPAR parameter values could be in comments
	else:
	    continue
  
	# Exclued parameter values which are in the UNUSED_PARAMS_TO_KEEP list
	try:
	    if inpa_mat in UNUSED_MAT_TO_KEEP:
		continue
	except KeyError:
	    pass
	
	# Grep for value in code. If the value doesn't appear than the input parameter might be unused
	# Check first part of files and if this failes check second part
	# If both fail than the argument doesn't exists  
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
	for ff in final_fail_printout:
	    print "The following inpar parameter only exists in drt_validmaterials.cpp"
	    owner = subprocess.check_output('grep ' + ff + ' ./praise.txt', shell=True)
	    owner = owner.strip(' )(\n,')
	    owner = owner[owner.index(' '):]
	    owner = (owner).strip()
	    owner = (owner.split(' '))
	    print ff, " ".join( ['' for i in range(55-len(ff))]),'last changed from: ', owner[0]
	sys.exit(1)	