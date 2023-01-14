set(OBJS_BINSTRATEGY
    # cmake-format: sortable
    binstrategy.cpp
    binstrategy_meshfree_bin.cpp
    binstrategy_meshfree_multibin.cpp
    binstrategy_utils.cpp
    )

baci_add_library(lib_binstrategy ${OBJS_BINSTRATEGY})

set(LIBS_LIB_BINSTRATEGY
    # cmake-format: sortable
    comm
    ele_bele
    ele_fluid
    ele_rigidsphere
    ele_scatra
    linalg
    rebalance
    )

target_link_libraries(lib_binstrategy PUBLIC ${LIBS_LIB_BINSTRATEGY})
