set(OBJS_DRT_LIB
    # cmake-format: sortable
    connectivity.cpp
    drt_assemblestrategy.cpp
    drt_condition.cpp
    drt_condition_selector.cpp
    drt_condition_utils.cpp
    drt_conditiondefinition.cpp
    drt_container.cpp
    drt_cubic_spline_interpolation.cpp
    drt_dirichletextractor.cpp
    drt_discret.cpp
    drt_discret_conditions.cpp
    drt_discret_evaluate.cpp
    drt_discret_faces.cpp
    drt_discret_fillcomplete.cpp
    drt_discret_hdg.cpp
    drt_discret_partition.cpp
    drt_discret_utils.cpp
    drt_discret_xfem.cpp
    drt_discret_xwall.cpp
    drt_dofset.cpp
    drt_dofset_base.cpp
    drt_dofset_definedmapping_wrapper.cpp
    drt_dofset_gidbased_wrapper.cpp
    drt_dofset_independent.cpp
    drt_dofset_independent_pbc.cpp
    drt_dofset_merged_wrapper.cpp
    drt_dofset_pbc.cpp
    drt_dofset_proxy.cpp
    drt_dofset_transparent.cpp
    drt_dofset_transparent_independent.cpp
    drt_element.cpp
    drt_elementdefinition.cpp
    drt_elementtype.cpp
    drt_exporter.cpp
    drt_function.cpp
    drt_function_library.cpp
    drt_function_manager.cpp
    drt_functionvariables.cpp
    drt_globalproblem.cpp
    drt_locsys.cpp
    drt_matchingoctree.cpp
    drt_node.cpp
    drt_parobject.cpp
    drt_parobjectfactory.cpp
    drt_parobjectregister.cpp
    drt_parser.cpp
    drt_periodicbc.cpp
    drt_resulttest.cpp
    drt_utils.cpp
    drt_utils_createdis.cpp
    drt_utils_discret.cpp
    drt_utils_discret_dbc.cpp
    drt_utils_discret_evaluate.cpp
    drt_utils_factory.cpp
    drt_utils_gid_vector.cpp
    drt_utils_materials.cpp
    drt_utils_parallel.cpp
    epetra_utils.cpp
    function_of_time.cpp
    tensor_transformation.cpp
    voigt_notation.cpp
    )

baci_add_library(drt_lib ${OBJS_DRT_LIB})

set(LIBS_DRT_LIB
    # cmake-format: sortable
    alg_particle_engine
    alg_poromultiphase_scatra
    alg_statinvanalysis
    drt_fem_general
    drt_lib_binstrategy
    drt_lib_geometric_search
    drt_lib_immersed
    drt_nurbs_discret
    ele_discsh3
    ele_elemag
    ele_to3
    headers
    linalg
    )

target_link_libraries(drt_lib PUBLIC ${LIBS_DRT_LIB})

# TODO this module requires a lot more includes that are currently included transitively from linalg

set(OBJS_DRT_IO_LIB
    # cmake-format: sortable
    drt_domainreader.cpp
    drt_elementreader.cpp
    drt_inputreader.cpp
    drt_linedefinition.cpp
    drt_materialdefinition.cpp
    drt_nodereader.cpp
    drt_particlereader.cpp
    )

baci_add_library(drt_lib_io ${OBJS_DRT_IO_LIB})

set(LIBS_DRT_LIB_IO
    # cmake-format: sortable
    drt_lib
    drt_lib_fiber
    )

target_link_libraries(drt_lib_io PUBLIC ${LIBS_DRT_LIB_IO})

add_subdirectory(error)
