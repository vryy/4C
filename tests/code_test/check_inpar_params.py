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

import sys, subprocess, re

# Dictionary of unused parameters we want to keep in the code
UNUSED_PARAMS_TO_KEEP = [ 
			# IrmLike isn't supported any longer since commit 17668
			'INPAR::STR::midavg_imrlike',		
			'INPAR::THR::midavg_imrlike', 
			# value if no solver is given
			'INPAR::STR::midavg_vague',
			'INPAR::STR::pred_vague',
			'INPAR::STR::soltech_vague',
			'INPAR::THR::midavg_vague',		
			'INPAR::THR::pred_vague',			
			'INPAR::THR::soltech_vague',
			# Parameters for path-following techniques, which currently aren't supported in Baci
			'INPAR::STR::control_arc1',
			'INPAR::STR::control_arc2',
			'INPAR::STR::control_disp',
			'INPAR::STR::control_load',
			# Paramters are actually used but could be replaced by a true/false flag
			'INPAR::SCATRA::evalmat_element_center',
			'INPAR::SCATRA::evaltau_element_center',
			'INPAR::SCATRA::penalty_method_none',
			'INPAR::ELCH::elch_mov_bndry_fully_transient',
                        # Parameter is necessary, but only used in validparameters.cpp
                        # NOTE To be implemented in future
                        'INPAR::CONTACT::wear_both_ale',
                        'INPAR::CONTACT::wear_both_map',
                        'INPAR::CONTACT::wear_archard', 
                        'INPAR::FLUID::EOS_CONV_CROSS_xfem_gp',
                        'INPAR::SOLVER::PA_AMG',                        
                        'INPAR::STR::stat_inv_mp_none',
                        'INPAR::XFEM::MSH_L2_Proj_part',
                        'INPAR::XFEM::interface_vel_init_zero',
                        'INPAR::CONTACT::bsm_cpp',
                        'INPAR::CONTACT::bsm_partially',
                        'INPAR::CONTACT::bsm_smoothed'
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
    
    subprocess.call('svn praise ' + global_src_path +  'drt_inpar/drt_validparameters.cpp > praise.txt' , shell=True )

    for inpa_para in progress('Searching inpar params', inpar_parameters):
    
	# Exclude commentary parts
	if not inpa_para.strip(' ')[:2] == '//':
    
	    # Get rid of leading and following brackets and of semicolons
	    inpa_para = re.sub(r'\(',',', inpa_para)
	    inpa_para = re.sub(r'\)',',', inpa_para)
	    inpa_para = re.sub(r'\;',',', inpa_para)
	    
	    # Search for keyword INPAR, if not found continue with loop
	    result_search = re.search('INPAR', inpa_para)
	    if result_search:
		inpa_para = inpa_para[result_search.start():]
		inpa_para = inpa_para.split(',')[0]
	    else:
	        continue
	
	# else path necessary, since the INPAR parameter values could be in comments
	else:
	    continue
  
	# Exclued parameter values which are in the UNUSED_PARAMS_TO_KEEP list
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
	for ff in final_fail_printout:
	    owner = subprocess.check_output('grep ' + ff + ' ./praise.txt', shell=True)
	    owner = owner.strip(' )(\n,')
	    owner = owner[owner.index(' '):]
	    owner = (owner).strip()
	    owner = (owner.split(' '))
	    print ff, " ".join( ['' for i in range(55-len(ff))]),'last changed from: ', owner[0]
	sys.exit(1)	
