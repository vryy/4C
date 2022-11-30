set(OBJS_DRT_FSI_XFEM
    fsi_xfem_fluid.cpp
    fsi_xfem_monolithic.cpp
    fsi_xfem_algorithm.cpp
    coupling_comm_manager.cpp
    XFScoupling_manager.cpp
    XFAcoupling_manager.cpp
    XFPcoupling_manager.cpp
    XFFcoupling_manager.cpp
    )

baci_add_library(alg_fsi_xfem ${OBJS_DRT_FSI_XFEM})
target_link_libraries(alg_fsi_xfem PUBLIC drt_adapter)
