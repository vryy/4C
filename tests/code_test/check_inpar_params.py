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

import sys, subprocess

# Dictionary of unused parameters we want to keep in the code
UNUSED_PARAMS_TO_KEEP = []


if __name__=='__main__':
  
    if len(sys.argv) < 2:
      print "usage: %s inputfile.dat source-path" % sys.argv[0]
      sys.exit(1)    
    
    # collect source files of src
    files_to_search = []
    global_src_path = sys.argv[1] + '/' + 'src/'
    
    source_headers = subprocess.check_output('ls ' + global_src_path, shell=True)
    for sh in source_headers.split():
	baci_files = subprocess.check_output('ls ' + global_src_path + sh, shell=True)
	baci_files = baci_files.split()
	files_to_search.extend( [global_src_path + sh + '/' + (baci_files)[i] for i in range(len(baci_files)) ] )
    
    # valid parameters file will be neclected
    files_to_search.remove(global_src_path + 'drt_inpar/drt_validparameters.cpp')
    
    print 'Start to search inpar params'
    
    # Initialization of Set for unused inpar parameters. Set used to avoid double entries
    fail = Set()
		
    inpar_parameters = subprocess.check_output('grep INPAR ' + global_src_path +  'drt_inpar/drt_validparameters.cpp', shell=True )
    inpar_parameters = inpar_parameters.split()
    
    for inpa_para in progress('Searching inpar params', inpar_parameters):
	
	inpa_para = inpa_para.strip(' \ntuple<int>(),;')
	inpa_para = (inpa_para.split(','))[0]
	
	if inpa_para[:5] == 'INPAR':
	  
	    try:
		if inpa_para in UNUSED_PARAMS_TO_KEEP:
		    continue
	    except KeyError:
		pass
	      
	    try:
		test = subprocess.check_output( '/bin/grep ' +  inpa_para + " " + " ".join(files_to_search), shell=True)
	    except subprocess.CalledProcessError: 
		fail.update([inpa_para])	    

    if not fail:
	print "Found no unused parameter"		
    else:
	# Retype to list in order to be able to sort
	final_fail_printout = [f for f in fail]
	final_fail_printout.sort()
	print "The following inpar parameter aren't used"
	print "\n".join(final_fail_printout)
	sys.exit(1)	