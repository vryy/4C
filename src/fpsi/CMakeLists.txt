set(OBJS_FPSI
    # cmake-format: sortable
    fpsi.cpp
    fpsi_coupling.cpp
    fpsi_dyn.cpp
    fpsi_monolithic.cpp
    fpsi_monolithic_plain.cpp
    fpsi_utils.cpp
    )

baci_add_library(alg_fpsi ${OBJS_FPSI})

set(LIBS_ALG_FPSI
    # cmake-format: sortable
    adapter
    )

baci_add_dependency(alg_fpsi ${LIBS_ALG_FPSI})
