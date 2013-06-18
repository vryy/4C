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
			  # Output parameters for GID will be kept. Output for SE (Strain Energy) will be implemented soon?!
			  'IO':                 ['OUTPUT_GID', 'STRUCT_SE'], 
			  # Path following paramters, which will be kept 
			  'STRUCTURAL DYNAMIC': ['CONTROLTYPE', 'CONTROLNODE']
			}


if __name__=='__main__':
  
    if len(sys.argv) != 3:
      print "usage: %s baci-exe source-path" % sys.argv[0]
      sys.exit(1)    
    
    # collect source files of src
    files_to_search = []
    name = sys.argv[2]
    
    if name[-1] == '/':
	global_src_path = sys.argv[2] + 'src/'
    else:	
	global_src_path = sys.argv[2] + '/' + 'src/'
    
    source_headers = subprocess.check_output('ls --hide=*.a ' + global_src_path, shell=True)
    for sh in source_headers.split():
	baci_files = subprocess.check_output('ls ' + global_src_path + sh, shell=True)
	baci_files = baci_files.split()
	files_to_search.extend( [global_src_path + sh + '/' + (baci_files)[i] for i in range(len(baci_files)) ] )
    
    # file drt_validparameters.cpp will be neclected in search
    files_to_search.remove(global_src_path + 'drt_inpar/drt_validparameters.cpp')
    
    # Create current default header file via the -d option of baci
    default_header_file = subprocess.check_output(sys.argv[1] + ' -d', shell=True)

    f = file('xxx_default.head','w')
    print >> f, default_header_file
    f.close()
    
    # Read in of default header
    section_names, sections = read_ccarat('xxx_default.head')
    
    print 'Start to search params'
    
    # fail dictionary
    fail = {}
    
    for section, section_values in sections.iteritems():
      
	if section == 'header':
	    continue

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
		
		# partitioning necessary due to maximal size of input arguments of bash console
		files_to_search_part1 = files_to_search[:len(files_to_search)/2 + 1]
		files_to_search_part2 = files_to_search[len(files_to_search)/2 + 1:]
		
		# Grep for value in code. If the value doesn't appear than the input parameter might be unused
		# Check first part of files and if this failes check second part
		# If both fail than the argument doesn't exists
		try:
		    test = subprocess.check_output( '/bin/grep ' +  sec_val[0] + " " + " ".join( files_to_search_part1 ), shell=True)
		except subprocess.CalledProcessError: 
		    try:
			test = subprocess.check_output( '/bin/grep ' +  sec_val[0] + " " + " ".join( files_to_search_part2 ), shell=True)
		    except subprocess.CalledProcessError: 
			if fail.has_key(section):
			    fail[section].update([sec_val[0]])
			else:
			    fail[section] = Set([sec_val[0]])
		
    if not fail:
	print "Found no unused input parameter"		
    else:
	for failsection, failvalues in fail.iteritems():
	    print 'Section', failsection, 'has the following input parameter, which only exists in drt_validparameters.cpp'
	    print "\n".join(failvalues)
	    print "\n"
	    
	sys.exit(1)    