#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Maintainer: A. Nagler
#
# Script for searching parameters which aren't used any longer

from read_ccarat_NIGHTLYTESTCASES import read_ccarat, write_ccarat
from elements    import bcdictionary, surfaces
from progress    import progress
from sets        import Set

import sys, subprocess, re

ELEMENTS_NOT_TO_SEARCH = []
MATERIALS_NOT_TO_SEARCH = []


if __name__=='__main__':
  
    if len(sys.argv) < 3:
      print "usage: %s baci-exe source-path" % sys.argv[0]
      sys.exit(1)    
    
    name = sys.argv[2]
    if name[-1] == '/':
	f = file(name + 'TestingFramework.cmake', 'r')
    else:
	f = file(name + '/' + 'TestingFramework.cmake', 'r')
    
    # Initialization of Set to avoid double counting
    nightly_testcases = Set()
    for line in f:
	if line[:9] == 'baci_test' and line[:20] != 'baci_test_Nested_Par':
	  
	    # remove baci_test
	    cur_line = line[9:]
	    # remove brackets
	    cur_line = cur_line.strip('()')
	    cur_line = cur_line.split(' ')
	    nightly_testcases.update( [cur_line[0] + '.dat'] )
	    
    f.close()
    
    # Retyping to list
    nightly_testcases = [n for n in nightly_testcases]
    
    # Create default header file via baci executible
    default_header_file = subprocess.check_output(sys.argv[1] + ' -d', shell=True)
    
    f = file('default.head','w')
    print >> f, default_header_file
    f.close()

    # Read in of default header
    section_names, sections = read_ccarat('default.head')
    
    
    # collect source files of src
    global_INPUT_path = sys.argv[2] + 'Input/'
    
    files_to_search = [ global_INPUT_path + nightly_testcases[i] for i in range(len(nightly_testcases)) ]
   
    # Check which materials have a nightly_testcases
    fail_mat = []
    for line in progress('Check materials',sections['MATERIALS']):
	if len(line)>2:
	    try:
		material = line[ line.index('//MAT') + 2 ]
	    except ValueError:
		continue
	      
	    if material in MATERIALS_NOT_TO_SEARCH:
		continue
	
	    try:
		test = subprocess.check_output('grep ' + material + ' ' + " ".join(files_to_search), shell=True)
	    except subprocess.CalledProcessError:
		fail_mat.append( material )
    
    # Check which element types have a nightly_testcases
    fail_ele = []
    element_sections = {'STRUCTURE ELEMENTS'		: 'structure elements', 
			'FLUID ELEMENTS'		: 'fluid elements', 
			'TRANSPORT ELEMENTS'		: 'transport elements', 
			'ALE ELEMENTS'			: 'ale elements', 
			'THERMO ELEMENTS'		: 'thermo elements', 
			'ARTERY ELEMENTS'		: 'artery elements', 
			'REDUCED D AIRWAYS ELEMENTS'	: 'red airways elements'
			}
    for ele_sec in element_sections.keys():
	
	for line in progress(element_sections[ele_sec],sections[ele_sec]):
	  
	    try:
		element = line[2] + ' ' + line[3]
	    except IndexError:
		continue
	      
	    if element in  ELEMENTS_NOT_TO_SEARCH:
		continue
	
	    try:
		test = subprocess.check_output('grep ' + "'" + element + "'" + ' ' + " ".join(files_to_search), shell=True)
	    except subprocess.CalledProcessError:
		fail_ele.append( element )		
    
    
    if not fail_mat and not fail_ele:
	print 'All implemented elements and materials have a nightly testcases'
    else:
	if fail_mat:
	    print 'The following materials don\'t have a nightly testcase'
	    print "\n"
	    print "\n".join(fail_mat)
    
 	if fail_ele:
	    print 'The following elements don\'t have a nightly testcase'
	    print "\n"
	    print "\n".join(fail_ele)   
	    
	sys.exit(1)    