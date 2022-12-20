set(OBJS_DISCSH3_LIB
    # cmake-format: sortable
    discsh3.cpp
    discsh3_evaluate.cpp
    discsh3_input.cpp
    discsh3_line.cpp
    discsh3_line_evaluate.cpp
    )

baci_add_library(ele_discsh3 ${OBJS_DISCSH3_LIB})

set(LIBS_ELE_DISCSH3
    # cmake-format: sortable
    linalg
    )

target_link_libraries(ele_discsh3 PUBLIC ${LIBS_ELE_DISCSH3})
