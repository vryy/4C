#
### General Setup

cmake_minimum_required(VERSION 2.8)
project(contact_module CXX C Fortran)

set(OBJS_DRT_LIB_MODULE
  src/drt_lib/drt_assemblestrategy.cpp
  src/drt_lib/drt_container.cpp
  src/drt_lib/drt_condition.cpp
  src/drt_lib/drt_condition_selector.cpp
  src/drt_lib/drt_condition_utils.cpp
  src/drt_lib/drt_discret.cpp
  src/drt_lib/drt_discret_fillcomplete.cpp
  src/drt_lib/drt_discret_partition.cpp
  src/drt_lib/drt_discret_evaluate.cpp
  src/drt_lib/drt_discret_conditions.cpp
  src/drt_lib/drt_discret_utils.cpp
  src/drt_lib/drt_dofset.cpp
  src/drt_lib/drt_dofset_base.cpp
  src/drt_lib/drt_dofset_proxy.cpp
  src/drt_lib/drt_dofset_subproxy.cpp
  src/drt_lib/drt_dofset_pbc.cpp
  src/drt_lib/drt_element.cpp
  src/drt_lib/drt_node.cpp
  src/drt_lib/drt_parobject.cpp
  src/drt_lib/drt_parobjectfactory.cpp
  src/drt_lib/drt_utils.cpp
  src/drt_lib/drt_utils_factory.cpp
  src/drt_lib/drt_utils_parmetis.cpp
  src/drt_lib/drt_utils_parallel.cpp
  src/drt_io/io_pstream.cpp
  src/drt_geometry/position_array.cpp
)
set(OBJS_DRT_CONTACT_LIB_MODULE
  src/drt_contact/contact_element.cpp
  src/drt_contact/contact_node.cpp
  src/drt_contact/friction_node.cpp
)
set(OBJS_DRT_UTILS_MODULE
  src/drt_fem_general/drt_utils_integration.cpp
  src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp
  src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp
  src/drt_fem_general/drt_utils_bspline.cpp
  src/drt_fem_general/drt_utils_boundary_integration.cpp
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_lib/drt_dserror.cpp
  src/drt_lib/drt_elementtype.cpp
  src/drt_lib/drt_exporter.cpp
)
set(OBJS_DRT_MORTAR_MODULE
  src/drt_mortar/mortar_binarytree.cpp
  src/drt_mortar/mortar_coupling2d.cpp
  src/drt_mortar/mortar_coupling3d.cpp
  src/drt_mortar/mortar_coupling3d_classes.cpp
  src/drt_mortar/mortar_utils.cpp
  src/drt_mortar/mortar_dofset.cpp
  src/drt_mortar/mortar_integrator.cpp
  src/drt_mortar/mortar_interface.cpp
  src/drt_mortar/mortar_projector.cpp
  src/drt_mortar/mortar_strategy_base.cpp
)
set(OBJS_DRT_MORTAR_LIB_MODULE
  src/drt_mortar/mortar_element.cpp
  src/drt_mortar/mortar_element_integrator.cpp
  src/drt_mortar/mortar_element_shapefct.cpp
  src/drt_mortar/mortar_node.cpp
)
set(OBJS_LINALG_MODULE
  src/linalg/linalg_blocksparsematrix.cpp
  src/linalg/linalg_mapextractor.cpp
  src/linalg/linalg_serialdensematrix.cpp
  src/linalg/linalg_serialdensevector.cpp
  src/linalg/linalg_sparsematrix.cpp
  src/linalg/linalg_sparsematrixbase.cpp
  src/linalg/linalg_utils.cpp
  src/linalg/linalg_multiply.cpp
)
set(OBJS_DRT_CONTACT_MODULE
  src/drt_contact/contact_abstract_strategy.cpp
  src/drt_contact/contact_coupling2d.cpp
  src/drt_contact/contact_coupling3d.cpp
  src/drt_contact/contact_interpolator.cpp
  src/drt_contact/contact_integrator.cpp
  src/drt_contact/contact_integrator_factory.cpp
  src/drt_contact/contact_integrator_assemble.cpp
  src/drt_contact/contact_interface.cpp
  src/drt_contact/contact_lagrange_strategy.cpp
  src/drt_contact/contact_penalty_strategy.cpp
  src/drt_contact/meshtying_abstract_strategy.cpp
  src/drt_contact/meshtying_lagrange_strategy.cpp
  src/drt_contact/meshtying_penalty_strategy.cpp
  src/drt_contact/selfcontact_binarytree.cpp
  src/drt_contact/contact_wear_interface.cpp
  src/drt_contact/contact_wear_lagrange_strategy.cpp
)
###------------------------------------------------------------------ Binaries
# Define the libraries. The link order is important!
add_library(alg_contact_mod ${OBJS_DRT_CONTACT_MODULE})
add_library(drt_mortar_mod ${OBJS_DRT_MORTAR_MODULE})

# libraries that are used
add_library(drt_lib_mod ${OBJS_DRT_LIB_MODULE})
add_library(drt_lib_contact_mod ${OBJS_DRT_CONTACT_LIB_MODULE} ${OBJS_DRT_BEAMCONTACT_LIB_MODULE})
add_library(drt_lib_mortar_mod ${OBJS_DRT_MORTAR_LIB_MODULE})
add_library(drt_utils_mod ${OBJS_DRT_UTILS_MODULE})
add_library(linalg_mod ${OBJS_LINALG_MODULE})

set(OBJS_TEST_MAIN
  src/main.cpp
)

# the test executable
add_executable (contact_test ${OBJS_TEST_MAIN})
set_target_properties (contact_test PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)

target_link_libraries (contact_test
alg_contact_mod
drt_mortar_mod
drt_lib_mortar_mod
drt_lib_contact_mod
drt_lib_mod
drt_utils_mod
linalg_mod
${LIBRARIES}
${LAPACK_LIBRARIES}
)

add_custom_command (
OUTPUT ${contact_module_SOURCE_DIR}/src/main.cpp
COMMAND ${contact_module_SOURCE_DIR}/generate_contact_module.py
)

add_custom_target(genContactModule ALL echo
DEPENDS ${contact_module_SOURCE_DIR}/src/main.cpp
)

add_dependencies(contact_test genContactModule)
