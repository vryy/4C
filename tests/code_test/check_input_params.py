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
    name = sys.argv[2]
    
    if name[-1] == '/':
	global_src_path = name + 'src/'
    else:	
	global_src_path = name + '/' + 'src/'
    
    # search for *.H and *.cpp files in the given src path	
    baci_heads = subprocess.check_output('find ' + global_src_path + ' -name *.H', shell=True)
    baci_heads = baci_heads.split()	
    
    baci_cpp   = subprocess.check_output('find ' + global_src_path + ' -name *.cpp', shell=True)
    baci_cpp   = baci_cpp.split()
    
    # exclude validparameters.H/.cpp from search
    baci_heads.remove(global_src_path + 'drt_inpar/drt_validparameters.H')
    baci_cpp.remove(global_src_path + 'drt_inpar/drt_validparameters.cpp')
    
    # Create current default header file via the -d option of baci
    default_header_file = subprocess.check_output(sys.argv[1] + ' -d', shell=True)

    f = file('xxx_default.head','w')
    print >> f, default_header_file
    f.close()
    
    # Read in of default header
    section_names, sections = read_ccarat('xxx_default.head')
    
    subprocess.call('svn praise ' + global_src_path +  'drt_inpar/drt_validparameters.cpp > praise.txt' , shell=True )    
    
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
		
		# Grep for value in code. If the value doesn't appear than the input parameter might be unused
		# Check first part of files and if this failes check second part
		# If both fail than the argument doesn't exists
		try:
		    test = subprocess.check_output( '/bin/grep ' +  sec_val[0] + " " + " ".join( baci_heads ), shell=True)
		except subprocess.CalledProcessError: 
		    try:
			test = subprocess.check_output( '/bin/grep ' +  sec_val[0] + " " + " ".join( baci_cpp ), shell=True)
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
	    for ff in failvalues:
		owner = subprocess.check_output('grep ' + ff + ' ./praise.txt', shell=True)
		owner = owner.strip(' )(\n,')
		owner = owner[owner.index(' '):]
		owner = (owner).strip()
		owner = (owner.split(' '))
		print ff, " ".join( ['' for i in range(55-len(ff))]),'last changed from: ', owner[0]		
    
	sys.exit(1)    