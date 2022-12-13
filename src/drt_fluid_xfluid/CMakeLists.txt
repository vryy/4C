set(OBJS_DRT_FLUID_XFLUID
    # cmake-format: sortable
    xfluid.cpp
    xfluid_functions.cpp
    xfluid_functions_combust.cpp
    xfluid_levelset_coupling_algorithm.cpp
    xfluid_outputservice.cpp
    xfluid_state.cpp
    xfluid_state_creator.cpp
    xfluidfluid.cpp
    xfluidfluid_state.cpp
    xfluidresulttest.cpp
    )

baci_add_library(alg_xfluid ${OBJS_DRT_FLUID_XFLUID})

set(LIBS_ALG_XFLUID
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(alg_xfluid PUBLIC ${LIBS_ALG_XFLUID})
