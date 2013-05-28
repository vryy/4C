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

# Dictionary of parameters which only appear in validparams
PARAMS_EXISTING_ONLY_IN_VALIDPARAMS = {
				      'TOPOLOGY OPTIMIZATION CONTROL': ['CONV_CHECK_TYPE']
				      }
# Dictionary of unused parameters we want to keep in the code
UNUSED_PARAMS_TO_KEEP = { 
			  'IO':                 ['OUTPUT_GID', 'STRUCT_SE'], 
			  'STRUCTURAL DYNAMIC': ['CONTROLTYPE', 'CONTROLNODE']
			}


if __name__=='__main__':
  
    if len(sys.argv) < 3:
      print "usage: %s inputfile.dat baci-exe source-path" % sys.argv[0]
      sys.exit(1)    
    
    # collect source files of src
    files_to_search = []
    global_src_path = sys.argv[2] + '/' + 'src/'
    
    source_headers = subprocess.check_output('ls ' + global_src_path, shell=True)
    for sh in source_headers.split():
	baci_files = subprocess.check_output('ls ' + global_src_path + sh, shell=True)
	baci_files = baci_files.split()
	files_to_search.extend( [global_src_path + sh + '/' + (baci_files)[i] for i in range(len(baci_files)) ] )
    
    # file drt_validparameters.cpp will be neclected in search
    files_to_search.remove(global_src_path + 'drt_inpar/drt_validparameters.cpp')
    
    # Create current default header file via the -d option of baci
    default_header_file = subprocess.check_output(sys.argv[1] + ' -d', shell=True)

    f = file('default.head','w')
    print >> f, default_header_file
    f.close()
    
    # Read in of default header
    section_names, sections = read_ccarat('default.head')
    
    print 'Start to search params'
    
    # fail dictionary
    fail = {}
    
    for section, section_values in sections.iteritems():

	for sec_val in section_values:
	  
	    # skip commentary parts of default.head
	    # NOTE case if empty lines have to be catched
	    if sec_val and sec_val[0][:2] != '//':
		    
		# skip values defined in dictionaries
		try:
		    if sec_val[0] in PARAMS_EXISTING_ONLY_IN_VALIDPARAMS[section]:
			continue
		except KeyError:
		    pass
		  
		try:
		    if sec_val[0] in UNUSED_PARAMS_TO_KEEP[section]:
			continue
		except KeyError:
		    pass		  
		
		# Grep for value in code. If the value doesn't appear than the input parameter might be unused
		try:
		    test = subprocess.check_output( '/bin/grep ' +  sec_val[0] + " " + " ".join(files_to_search), shell=True)
		except subprocess.CalledProcessError: 
		    #print 'ERROR: Input Parameter', sec_val[0], 'of section', section, 'only found in drt_validparameters.cpp' 
		    #fail = True
		    if fail.has_key(section):
			fail[section].update([sec_val[0]])
		    else:
			fail[section] = Set([sec_val[0]])
		
    if not fail:
	print "Found no unused parameter"		
    else:
	for failsection, failvalues in fail.iteritems():
	    print 'Section', failsection, 'has the following unused input parameter'
	    print "\n".join(failvalues)
	    print "\n"
	    
	sys.exit(1)    