set(OBJS_DRT_BINSTRATEGY
    binning_strategy.cpp binning_strategy_utils.cpp drt_meshfree_bin.cpp drt_meshfree_multibin.cpp
    )

baci_add_library(drt_lib_binstrategy ${OBJS_DRT_BINSTRATEGY})
target_link_libraries(drt_lib_binstrategy PUBLIC linalg)
