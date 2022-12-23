set(OBJS_FSI_XFEM
    # cmake-format: sortable
    coupling_comm_manager.cpp
    fsi_xfem_algorithm.cpp
    fsi_xfem_fluid.cpp
    fsi_xfem_monolithic.cpp
    XFAcoupling_manager.cpp
    XFFcoupling_manager.cpp
    XFPcoupling_manager.cpp
    XFScoupling_manager.cpp
    )

baci_add_library(alg_fsi_xfem ${OBJS_FSI_XFEM})

set(LIBS_ALG_FSI_XFEM
    # cmake-format: sortable
    adapter
    )

target_link_libraries(alg_fsi_xfem PUBLIC ${LIBS_ALG_FSI_XFEM})
