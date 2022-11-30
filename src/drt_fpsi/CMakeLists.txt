set(OBJS_DRT_FPSI
    fpsi_utils.cpp
    fpsi_dyn.cpp
    fpsi_monolithic.cpp
    fpsi_monolithic_plain.cpp
    fpsi.cpp
    fpsi_coupling.cpp
    )

baci_add_library(alg_fpsi ${OBJS_DRT_FPSI})
target_link_libraries(alg_fpsi PUBLIC drt_adapter)
