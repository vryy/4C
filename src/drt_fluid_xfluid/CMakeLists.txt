set(OBJS_DRT_FLUID_XFLUID
    xfluid.cpp
    xfluidfluid.cpp
    xfluidfluid_state.cpp
    xfluidresulttest.cpp
    xfluid_state.cpp
    xfluid_state_creator.cpp
    xfluid_levelset_coupling_algorithm.cpp
    xfluid_outputservice.cpp
    )

baci_add_library(alg_xfluid ${OBJS_DRT_FLUID_XFLUID})
target_link_libraries(alg_xfluid PUBLIC drt_lib)
