set(OBJS_DRT_XFEM
    xfem_edgestab.cpp
    xfem_discretization_utils.cpp
    xfem_dofset_transparent_independent.cpp
    xfem_dofset.cpp
    xfem_mesh_projector.cpp
    xfem_neumann.cpp
    xfem_interface_utils.cpp
    xfem_utils.cpp
    xfield_state.cpp
    xfield_field_coupling_dofset.cpp
    xfield_field_coupling.cpp
    xfem_multi_field_mapextractor.cpp
    xfluid_timeInt_base.cpp
    xfluid_timeInt_std_SemiLagrange.cpp
    xfluid_timeInt.cpp
    )

baci_add_library(drt_lib_xfem ${OBJS_DRT_XFEM})
target_link_libraries(drt_lib_xfem PUBLIC linalg drt_lib drt_inpar)

set(OBJS_DRT_XFEM_MANAG_LIB
    xfem_coupling_levelset.cpp
    xfem_coupling_base.cpp
    xfem_coupling_mesh.cpp
    xfem_coupling_fpi_mesh.cpp
    xfem_coupling_mesh_coupled_levelset.cpp
    xfem_condition_manager.cpp
    xfem_xfluid_contact_communicator.cpp
    )

baci_add_library(drt_lib_xfem_manag ${OBJS_DRT_XFEM_MANAG_LIB})
target_link_libraries(drt_lib_xfem_manag PUBLIC linalg drt_lib drt_inpar)
