set(OBJS_DRT_XFEM
    # cmake-format: sortable
    xfem_discretization_utils.cpp
    xfem_dofset.cpp
    xfem_dofset_transparent_independent.cpp
    xfem_edgestab.cpp
    xfem_interface_utils.cpp
    xfem_mesh_projector.cpp
    xfem_multi_field_mapextractor.cpp
    xfem_neumann.cpp
    xfem_utils.cpp
    xfield_field_coupling.cpp
    xfield_field_coupling_dofset.cpp
    xfield_state.cpp
    xfluid_timeInt.cpp
    xfluid_timeInt_base.cpp
    xfluid_timeInt_std_SemiLagrange.cpp
    )

baci_add_library(drt_lib_xfem ${OBJS_DRT_XFEM})

set(LIBS_DRT_LIB_XFEM
    # cmake-format: sortable
    drt_inpar
    drt_lib
    linalg
    )

target_link_libraries(drt_lib_xfem PUBLIC ${LIBS_DRT_LIB_XFEM})

set(OBJS_DRT_XFEM_MANAG_LIB
    # cmake-format: sortable
    xfem_condition_manager.cpp
    xfem_coupling_base.cpp
    xfem_coupling_fpi_mesh.cpp
    xfem_coupling_levelset.cpp
    xfem_coupling_mesh.cpp
    xfem_coupling_mesh_coupled_levelset.cpp
    xfem_xfluid_contact_communicator.cpp
    )

baci_add_library(drt_lib_xfem_manag ${OBJS_DRT_XFEM_MANAG_LIB})

set(LIBS_DRT_LIB_XFEM_MANAG
    # cmake-format: sortable
    drt_inpar
    drt_lib
    linalg
    )

target_link_libraries(drt_lib_xfem_manag PUBLIC ${LIBS_DRT_LIB_XFEM_MANAG})
