set(OBJS_FLUID_XFLUID
    # cmake-format: sortable
    fluid_xfluid.cpp
    fluid_xfluid_fluid.cpp
    fluid_xfluid_fluid_state.cpp
    fluid_xfluid_functions.cpp
    fluid_xfluid_functions_combust.cpp
    fluid_xfluid_levelset_coupling_algorithm.cpp
    fluid_xfluid_outputservice.cpp
    fluid_xfluid_resulttest.cpp
    fluid_xfluid_state.cpp
    fluid_xfluid_state_creator.cpp
    )

baci_add_library(alg_xfluid ${OBJS_FLUID_XFLUID})

set(LIBS_ALG_XFLUID
    # cmake-format: sortable
    lib
    )

target_link_libraries(alg_xfluid PUBLIC ${LIBS_ALG_XFLUID})
