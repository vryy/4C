set(OBJS_FS3I
    # cmake-format: sortable
    ac_fsi.cpp
    ac_fsi_lts.cpp
    biofilm_fsi.cpp
    biofilm_fsi_utils.cpp
    fps3i_partitioned.cpp
    fps3i_partitioned_1wc.cpp
    fs3i.cpp
    fs3i_dyn.cpp
    fs3i_partitioned.cpp
    fs3i_partitioned_1wc.cpp
    fs3i_partitioned_2wc.cpp
    )

baci_add_library(alg_fs3i ${OBJS_FS3I})

set(LIBS_ALG_FS3I
    # cmake-format: sortable
    inpar
    )

target_link_libraries(alg_fs3i PUBLIC ${LIBS_ALG_FS3I})
