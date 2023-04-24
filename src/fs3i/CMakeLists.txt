set(OBJS_FS3I
    # cmake-format: sortable
    fs3i.cpp
    fs3i_ac_fsi.cpp
    fs3i_ac_fsi_lts.cpp
    fs3i_biofilm_fsi.cpp
    fs3i_biofilm_fsi_utils.cpp
    fs3i_dyn.cpp
    fs3i_fps3i_partitioned.cpp
    fs3i_fps3i_partitioned_1wc.cpp
    fs3i_partitioned.cpp
    fs3i_partitioned_1wc.cpp
    fs3i_partitioned_2wc.cpp
    )

baci_add_library(alg_fs3i ${OBJS_FS3I})

set(LIBS_ALG_FS3I
    # cmake-format: sortable
    inpar
    )

baci_add_dependency(alg_fs3i ${LIBS_ALG_FS3I})
