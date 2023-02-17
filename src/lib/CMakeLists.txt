set(OBJS_LIB
    # cmake-format: sortable
    lib_assemblestrategy.cpp
    lib_condition.cpp
    lib_condition_selector.cpp
    lib_condition_utils.cpp
    lib_conditiondefinition.cpp
    lib_connectivity.cpp
    lib_container.cpp
    lib_cubic_spline_interpolation.cpp
    lib_dirichletextractor.cpp
    lib_discret.cpp
    lib_discret_broadcast.cpp
    lib_discret_conditions.cpp
    lib_discret_evaluate.cpp
    lib_discret_faces.cpp
    lib_discret_fillcomplete.cpp
    lib_discret_hdg.cpp
    lib_discret_partition.cpp
    lib_discret_utils.cpp
    lib_discret_xfem.cpp
    lib_discret_xwall.cpp
    lib_dofset.cpp
    lib_dofset_base.cpp
    lib_dofset_definedmapping_wrapper.cpp
    lib_dofset_gidbased_wrapper.cpp
    lib_dofset_independent.cpp
    lib_dofset_independent_pbc.cpp
    lib_dofset_merged_wrapper.cpp
    lib_dofset_pbc.cpp
    lib_dofset_proxy.cpp
    lib_dofset_transparent.cpp
    lib_dofset_transparent_independent.cpp
    lib_element.cpp
    lib_elementdefinition.cpp
    lib_elementtype.cpp
    lib_epetra_utils.cpp
    lib_exporter.cpp
    lib_function.cpp
    lib_function_library.cpp
    lib_function_manager.cpp
    lib_function_of_time.cpp
    lib_functionvariables.cpp
    lib_globalproblem.cpp
    lib_gridgenerator.cpp
    lib_locsys.cpp
    lib_matchingoctree.cpp
    lib_node.cpp
    lib_parobject.cpp
    lib_parobjectfactory.cpp
    lib_parobjectregister.cpp
    lib_periodicbc.cpp
    lib_resulttest.cpp
    lib_symbolic_expression.cpp
    lib_tensor_transformation.cpp
    lib_utils.cpp
    lib_utils_createdis.cpp
    lib_utils_discret.cpp
    lib_utils_discret_dbc.cpp
    lib_utils_discret_evaluate.cpp
    lib_utils_factory.cpp
    lib_utils_gid_vector.cpp
    lib_utils_materials.cpp
    lib_utils_parallel.cpp
    lib_voigt_notation.cpp
    )

baci_add_library(lib ${OBJS_LIB})

set(LIBS_LIB
    # cmake-format: sortable
    alg_particle_engine
    alg_poromultiphase_scatra
    alg_statinvanalysis
    ele_discsh3
    ele_elemag
    ele_to3
    fem_general
    headers
    lib_binstrategy
    lib_geometric_search
    lib_immersed
    linalg
    nurbs_discret
    )

baci_add_dependency(lib ${LIBS_LIB})

# TODO this module requires a lot more includes that are currently included transitively from linalg

set(OBJS_IO_LIB
    # cmake-format: sortable
    lib_domainreader.cpp
    lib_elementreader.cpp
    lib_inputreader.cpp
    lib_linedefinition.cpp
    lib_materialdefinition.cpp
    lib_meshreader.cpp
    lib_nodereader.cpp
    lib_particlereader.cpp
    )

baci_add_library(lib_io ${OBJS_IO_LIB})

set(LIBS_LIB_IO
    # cmake-format: sortable
    lib
    lib_fiber
    )

baci_add_dependency(lib_io ${LIBS_LIB_IO})

add_subdirectory(error)
