set(OBJS_DRT_LIB
    drt_periodicbc.cpp
    drt_assemblestrategy.cpp
    drt_condition.cpp
    drt_conditiondefinition.cpp
    drt_container.cpp
    drt_cubic_spline_interpolation.cpp
    drt_discret.cpp
    drt_discret_faces.cpp
    drt_discret_hdg.cpp
    drt_discret_xwall.cpp
    drt_discret_xfem.cpp
    drt_discret_conditions.cpp
    drt_discret_evaluate.cpp
    drt_discret_fillcomplete.cpp
    drt_discret_partition.cpp
    drt_discret_utils.cpp
    drt_dofset.cpp
    drt_dofset_base.cpp
    drt_dofset_proxy.cpp
    drt_dofset_gidbased_wrapper.cpp
    drt_dofset_merged_wrapper.cpp
    drt_dofset_definedmapping_wrapper.cpp
    drt_dofset_transparent.cpp
    drt_dofset_transparent_independent.cpp
    drt_dofset_pbc.cpp
    drt_dofset_independent.cpp
    drt_dofset_independent_pbc.cpp
    drt_element.cpp
    drt_function.cpp
    drt_functionvariables.cpp
    drt_function_library.cpp
    drt_function_manager.cpp
    function_of_time.cpp
    drt_globalproblem.cpp
    drt_locsys.cpp
    drt_node.cpp
    drt_parobject.cpp
    drt_parobjectfactory.cpp
    drt_utils_createdis.cpp
    drt_utils_materials.cpp
    drt_utils_factory.cpp
    drt_utils_gid_vector.cpp
    drt_utils_discret.cpp
    drt_utils_discret_dbc.cpp
    drt_utils_discret_evaluate.cpp
    connectivity.cpp
    tensor_transformation.cpp
    )

baci_add_library(drt_lib ${OBJS_DRT_LIB})
target_link_libraries(drt_lib PUBLIC linalg drt_lib_binstrategy drt_lib_geometric_search)
# TODO this module requires a lot more includes that are currently included transitively from linalg

add_subdirectory(error)
