set(OBJS_T3CL_LIB
    # cmake-format: sortable
    truss3cl.cpp
    truss3cl_evaluate.cpp
    truss3cl_input.cpp
    )

baci_add_library(ele_t3cl ${OBJS_T3CL_LIB})

set(LIBS_ELET3CL
    # cmake-format: sortable
    lib
    )

baci_add_dependency(ele_t3cl ${LIBS_ELET3CL})
