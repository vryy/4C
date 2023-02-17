set(OBJS_FSI_XFEM
    # cmake-format: sortable
    fsi_xfem_algorithm.cpp
    fsi_xfem_coupling_comm_manager.cpp
    fsi_xfem_fluid.cpp
    fsi_xfem_monolithic.cpp
    fsi_xfem_XFAcoupling_manager.cpp
    fsi_xfem_XFFcoupling_manager.cpp
    fsi_xfem_XFPcoupling_manager.cpp
    fsi_xfem_XFScoupling_manager.cpp
    )

baci_add_library(alg_fsi_xfem ${OBJS_FSI_XFEM})

set(LIBS_ALG_FSI_XFEM
    # cmake-format: sortable
    adapter
    )

baci_add_dependency(alg_fsi_xfem ${LIBS_ALG_FSI_XFEM})
