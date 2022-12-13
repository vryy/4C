set(OBJS_DRT_FPSI
    # cmake-format: sortable
    fpsi.cpp
    fpsi_coupling.cpp
    fpsi_dyn.cpp
    fpsi_monolithic.cpp
    fpsi_monolithic_plain.cpp
    fpsi_utils.cpp
    )

baci_add_library(alg_fpsi ${OBJS_DRT_FPSI})

set(LIBS_ALG_FPSI
    # cmake-format: sortable
    drt_adapter
    )

target_link_libraries(alg_fpsi PUBLIC ${LIBS_ALG_FPSI})
