set(OBJS_T3_LIB
    # cmake-format: sortable
    truss3.cpp
    truss3_evaluate.cpp
    truss3_input.cpp
    truss3_scatra.cpp
    )

baci_add_library(ele_t3 ${OBJS_T3_LIB})

set(LIBS_ELET3
    # cmake-format: sortable
    lib
    )

baci_add_dependency(ele_t3 ${LIBS_ELET3})
