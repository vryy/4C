set(OBJS_DRT_T3CL_LIB
    # cmake-format: sortable
    truss3cl.cpp
    truss3cl_evaluate.cpp
    truss3cl_input.cpp
    )

baci_add_library(ele_t3cl ${OBJS_DRT_T3CL_LIB})

set(LIBS_ELET3CL
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(ele_t3cl PUBLIC ${LIBS_ELET3CL})
