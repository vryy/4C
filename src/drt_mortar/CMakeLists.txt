set(OBJS_DRT_MORTAR
    # cmake-format: sortable
    base_binarytree.cpp
    mortar_binarytree.cpp
    mortar_coupling2d.cpp
    mortar_coupling3d.cpp
    mortar_coupling3d_classes.cpp
    mortar_dofset.cpp
    mortar_integrator.cpp
    mortar_interface.cpp
    mortar_interface_tools.cpp
    mortar_interface_utils.cpp
    mortar_manager_base.cpp
    mortar_matrix_transform.cpp
    mortar_multifield_coupling.cpp
    mortar_projector.cpp
    mortar_strategy_base.cpp
    mortar_strategy_factory.cpp
    mortar_utils.cpp
    )

baci_add_library(drt_mortar ${OBJS_DRT_MORTAR})

set(LIBS_DRT_LIB_MORTAR
    # cmake-format: sortable
    drt_lib
    nlnsolver_nox
    )

target_link_libraries(drt_mortar PUBLIC ${LIBS_DRT_LIB_MORTAR})

set(OBJS_DRT_LIB_MORTAR
    # cmake-format: sortable
    mortar_element.cpp
    mortar_element_integrator.cpp
    mortar_element_shapefct.cpp
    mortar_node.cpp
    )

baci_add_library(drt_lib_mortar ${OBJS_DRT_LIB_MORTAR})

set(LIBS_DRT_LIB_MORTAR
    # cmake-format: sortable
    drt_lib
    nlnsolver_nox
    )

target_link_libraries(drt_lib_mortar PUBLIC ${LIBS_DRT_LIB_MORTAR})
