#!/usr/bin/python

#import os.path
import os
import sys
import math
import getpass

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    cmd = "rm -rf "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def createDir(path):
    """deletes the path entirely"""
    cmd = "mkdir "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def runCommand(cmd):
    """deletes the path entirely"""
    result = getstatusoutput(cmd)
    #if(result[0]!=0):
    #    raise RuntimeError(result[1])
    return result[1]

###########
# MAIN routine
def main(argv=None):

  # This is the destination folder for the contact module
  folder = "./"
  srcfolder = folder + "src/" # This is the destination folder for the source files

  # This is the BACI input folder
  inputfolder = "../../"

  # clean old files
  deleteDir(srcfolder)
  createDir(srcfolder)

  # copy files
  createDir(srcfolder + "drt_contact")
  runCommand("cp -a " + inputfolder + "src/drt_contact/* " + srcfolder + "drt_contact/.")

  createDir(srcfolder + "drt_mortar")
  runCommand("cp -a " + inputfolder + "src/drt_mortar/* " + srcfolder + "drt_mortar/.")

  ## copy single files
  createDir(srcfolder + "drt_inpar")
  runCommand("cp -a " + inputfolder + "src/drt_inpar/inpar_mortar.H " + srcfolder + "drt_inpar/.")
  runCommand("cp -a " + inputfolder + "src/drt_inpar/inpar_parameterlist_utils.H " + srcfolder + "drt_inpar/.")
  runCommand("cp -a " + inputfolder + "src/drt_inpar/inpar_contact.H " + srcfolder + "drt_inpar/.")
  runCommand("cp -a " + inputfolder + "src/drt_inpar/inpar_wear.H " + srcfolder + "drt_inpar/.")

  createDir(srcfolder + "drt_io")
  runCommand("cp -a " + inputfolder + "src/drt_io/io_pstream.H " + srcfolder + "drt_io/.")
  runCommand("cp -a " + inputfolder + "src/drt_io/io_pstream.cpp " + srcfolder + "drt_io/.")

  createDir(srcfolder + "drt_lib")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dserror.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dserror.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_node.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_node.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_parobject.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_parobject.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_pack_buffer.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_container.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_container.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_parobjectfactory.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_parobjectfactory.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_singletondestruction.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_fillcomplete.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_partition.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_conditions.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_evaluate.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_iterator.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_discret_utils.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_assemblestrategy.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_assemblestrategy.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition_selector.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition_selector.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition_utils.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_condition_utils.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_element.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_element.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_base.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_base.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_proxy.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_proxy.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_subproxy.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_subproxy.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_aux_proxy.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_pbc.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_dofset_pbc.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_exporter.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_exporter.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_factory.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_factory.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_elementtype.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_globalproblem_enums.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_elementdefinition.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_elementtype.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_colors.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_element_integration_select.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_parmetis.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_parmetis.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_parallel.H " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_utils_parallel.cpp " + srcfolder + "drt_lib/.")
  runCommand("cp -a " + inputfolder + "src/drt_lib/drt_elements_paramsinterface.H " + srcfolder + "drt_lib/.")

  createDir(srcfolder + "drt_geometry")
  runCommand("cp -a " + inputfolder + "src/drt_geometry/position_array.H " + srcfolder + "drt_geometry/.")
  runCommand("cp -a " + inputfolder + "src/drt_geometry/position_array.cpp " + srcfolder + "drt_geometry/.")

  createDir(srcfolder + "drt_fem_general")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_fem_shapefunctions.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_local_connectivity_matrices.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_integration.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_integration.cpp " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_nurbs_shapefunctions.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_bspline.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_bspline.cpp " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_boundary_integration.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_boundary_integration.cpp " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_gausspoints.H " + srcfolder + "drt_fem_general/.")
  runCommand("cp -a " + inputfolder + "src/drt_fem_general/drt_utils_gausspoints.cpp " + srcfolder + "drt_fem_general/.")

  createDir(srcfolder + "linalg")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_serialdensematrix.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_serialdensematrix.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_serialdensevector.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_serialdensevector.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_fixedsizematrix.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_fixedsizematrix.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_sparsematrix.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_sparsematrix.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_sparsematrixbase.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_sparsematrixbase.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_sparseoperator.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_mapextractor.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_mapextractor.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_utils.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_utils.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_blocksparsematrix.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_blocksparsematrix.cpp " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_multiply.H " + srcfolder + "linalg/.")
  runCommand("cp -a " + inputfolder + "src/linalg/linalg_multiply.cpp " + srcfolder + "linalg/.")

  createDir(srcfolder + "headers")
  runCommand("cp -a " + inputfolder + "src/headers/pairedvector.H " + srcfolder + "headers/.")
  runCommand("cp -a " + inputfolder + "src/headers/compiler_definitions.h " + srcfolder + "headers/.")

  createDir(srcfolder + "drt_structure_new")
  runCommand("cp -a " + inputfolder + "src/drt_structure_new/str_enum_lists.H " + srcfolder + "drt_structure_new/.")

  createDir(srcfolder + "solver_nonlin_nox")
  runCommand("cp -a " + inputfolder + "src/solver_nonlin_nox/nox_nln_constraint_interface_preconditioner.H " + srcfolder + "solver_nonlin_nox/.")
  runCommand("cp -a " + inputfolder + "src/solver_nonlin_nox/nox_nln_enum_lists.H " + srcfolder + "solver_nonlin_nox/.")

  ## copy all patch files
  patchfolder = folder + "patches/"
  #runCommand("cp -a " + patchfolder + "*.patch " + folder)
  ##deleteDir(patchfolder)
  ##createDir(patchfolder)
  ##runCommand("cp -a " + inputfolder + "*.patch " + patchfolder)
  ##runCommand("cp -a " + inputfolder + "*.patch " + folder)

  # apply patch files
  runCommand("for i in ./patches/*.patch; do patch -p0 < $i; done") # only GIT patches
  runCommand("for i in ./patches/*.patch; do patch -p1 < $i; done") # only GIT patches

  # copy test program to main folder
  runCommand ("cp -a tmpl/main.cpp src/main.cpp")
  runCommand ("cp -a tmpl/validParameters.H src/validParameters.H")

if __name__ == "__main__":
  sys.exit(main())

