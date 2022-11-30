set(OBJS_DRT_FS3I
    fs3i.cpp
    fs3i_dyn.cpp
    fs3i_partitioned.cpp
    fps3i_partitioned.cpp
    fs3i_partitioned_1wc.cpp
    fps3i_partitioned_1wc.cpp
    fs3i_partitioned_2wc.cpp
    ac_fsi.cpp
    ac_fsi_lts.cpp
    biofilm_fsi.cpp
    biofilm_fsi_utils.cpp
    )

baci_add_library(alg_fs3i ${OBJS_DRT_FS3I})
target_link_libraries(alg_fs3i PUBLIC drt_lib_inpar)
