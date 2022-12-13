set(OBJS_DRT_BINSTRATEGY
    # cmake-format: sortable
    binning_strategy.cpp
    binning_strategy_utils.cpp
    drt_meshfree_bin.cpp
    drt_meshfree_multibin.cpp
    )

baci_add_library(drt_lib_binstrategy ${OBJS_DRT_BINSTRATEGY})

set(LIBS_DRT_LIB_BINSTRATEGY
    # cmake-format: sortable
    drt_comm
    ele_bele
    ele_fluid
    ele_rigidsphere
    ele_scatra
    linalg
    rebalance
    )

target_link_libraries(drt_lib_binstrategy PUBLIC ${LIBS_DRT_LIB_BINSTRATEGY})
