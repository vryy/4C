#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Maintainer: A. Nagler
#
# Script for searching parameters which aren't used any longer

from read_ccarat_NIGHTLYTESTCASES import read_ccarat, write_ccarat
from elements    		  import bcdictionary, surfaces
from progress    		  import progress
from sets        		  import Set
from copy 			  import copy

import sys, subprocess, re

# Dictionary of unused parameters we want to keep in the code
UNUSED_PARAMS_TO_KEEP = [ 
			# IrmLike isn't supported any longer since commit 17668
			'INPAR::STR::midavg_imrlik',		
			'INPAR::THR::midavg_imrlik', 
			# value if no solver is given
			'INPAR::STR::midavg_vag',
			'INPAR::STR::pred_vag',
			'INPAR::STR::soltech_vag',
			'INPAR::THR::midavg_vag',		
			'INPAR::THR::pred_vag',			
			'INPAR::THR::soltech_vag',
			# Parameters for path-following techniques, which currently aren't supported in Baci
			'INPAR::STR::control_arc1',
			'INPAR::STR::control_arc2',
			'INPAR::STR::control_dis',
			'INPAR::STR::control_load',
			# Paramters are actually used but could be replaced by a true/false flag
			'INPAR::SCATRA::evalmat_element_center',
			'INPAR::SCATRA::evaltau_element_center',
			'INPAR::SCATRA::penalty_method_no',
			'INPAR::ELCH::elch_mov_bndry_fully_trans',
                        # Cavitation parameters to be used in the future
                        'INPAR::CAVITATION::TwoWayMomentum',
                        'INPAR::CAVITATION::TwoWayFull'
			]

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
    
    # exclude validparameters.H/.cpp from search
    baci_heads.remove(global_src_path + 'drt_inpar/drt_validparameters.H')
    baci_cpp.remove(global_src_path + 'drt_inpar/drt_validparameters.cpp')
    
    print 'Start to search inpar params'
    
    # Initialization of Set for unused inpar parameters. Set used to avoid double entries
    fail = Set()
		
    inpar_parameters = subprocess.check_output('grep INPAR ' + global_src_path +  'drt_inpar/drt_validparameters.cpp', shell=True )
    inpar_parameters = inpar_parameters.split()
    
    for inpa_para in progress('Searching inpar params', inpar_parameters):

	# special treatment of lines which include tuple<int> before the parameter
	# i.e. tuple<int>(INPAR::CAVIATION::TwoWayFull
	# method is to substitute the brackets ( with an ,
	inpa_para = re.sub(r'\(',',', inpa_para)
	inpa_para = (inpa_para.split(','))
	
	# search for desired inpar value
	for ip in inpa_para:
	    ip = ip.strip(' \n,;()')
	    if ip[:5]== 'INPAR':
		inpa_para = ip
		break
	
	# else path necessary, since the INPAR parameter values could be in comments
	else:
	    continue
  
	try:
	    if inpa_para in UNUSED_PARAMS_TO_KEEP:
		continue
	except KeyError:
	    pass
	
	# Grep for value in code. If the value doesn't appear than the input parameter might be unused
	# Check first part of files and if this failes check second part
	# If both fail than the argument doesn't exists  
	try:
	    test = subprocess.check_output( '/bin/grep ' +  inpa_para + " " + " ".join(baci_heads), shell=True)
	except subprocess.CalledProcessError: 
	    try:
		test = subprocess.check_output( '/bin/grep ' +  inpa_para + " " + " ".join(baci_cpp), shell=True)
	    except subprocess.CalledProcessError: 
		fail.update([inpa_para])	
		
    if not fail:
	print "Found no unused parameter"		
    else:
	# Retype to list in order to be able to sort
	final_fail_printout = [f for f in fail]
	final_fail_printout.sort()
	print "The following inpar parameter only exists in drt_validparameters.cpp"
	print "\n".join(final_fail_printout)
	sys.exit(1)	
