set(OBJS_LIB
    # cmake-format: sortable
    assemblestrategy.cpp
    condition.cpp
    condition_selector.cpp
    condition_utils.cpp
    conditiondefinition.cpp
    connectivity.cpp
    container.cpp
    cubic_spline_interpolation.cpp
    dirichletextractor.cpp
    discret.cpp
    discret_broadcast.cpp
    discret_conditions.cpp
    discret_evaluate.cpp
    discret_faces.cpp
    discret_fillcomplete.cpp
    discret_hdg.cpp
    discret_partition.cpp
    discret_utils.cpp
    discret_xfem.cpp
    discret_xwall.cpp
    dofset.cpp
    dofset_base.cpp
    dofset_definedmapping_wrapper.cpp
    dofset_gidbased_wrapper.cpp
    dofset_independent.cpp
    dofset_independent_pbc.cpp
    dofset_merged_wrapper.cpp
    dofset_pbc.cpp
    dofset_proxy.cpp
    dofset_transparent.cpp
    dofset_transparent_independent.cpp
    element.cpp
    elementdefinition.cpp
    elementtype.cpp
    epetra_utils.cpp
    exporter.cpp
    function.cpp
    function_library.cpp
    function_manager.cpp
    function_of_time.cpp
    functionvariables.cpp
    globalproblem.cpp
    locsys.cpp
    matchingoctree.cpp
    node.cpp
    parobject.cpp
    parobjectfactory.cpp
    parobjectregister.cpp
    periodicbc.cpp
    resulttest.cpp
    symbolic_expression.cpp
    tensor_transformation.cpp
    utils.cpp
    utils_createdis.cpp
    utils_discret.cpp
    utils_discret_dbc.cpp
    utils_discret_evaluate.cpp
    utils_factory.cpp
    utils_gid_vector.cpp
    utils_materials.cpp
    utils_parallel.cpp
    voigt_notation.cpp
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

target_link_libraries(lib PUBLIC ${LIBS_LIB})

# TODO this module requires a lot more includes that are currently included transitively from linalg

set(OBJS_IO_LIB
    # cmake-format: sortable
    domainreader.cpp
    elementreader.cpp
    inputreader.cpp
    linedefinition.cpp
    materialdefinition.cpp
    nodereader.cpp
    particlereader.cpp
    )

baci_add_library(lib_io ${OBJS_IO_LIB})

set(LIBS_LIB_IO
    # cmake-format: sortable
    lib
    lib_fiber
    )

target_link_libraries(lib_io PUBLIC ${LIBS_LIB_IO})

add_subdirectory(error)
