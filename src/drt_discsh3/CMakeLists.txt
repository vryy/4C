set(OBJS_DRT_DISCSH3_LIB
    discsh3.cpp
    discsh3_line.cpp
    discsh3_evaluate.cpp
    discsh3_line_evaluate.cpp
    discsh3_input.cpp
    )

baci_add_library(ele_discsh3 ${OBJS_DRT_DISCSH3_LIB})

target_link_libraries(ele_discsh3 PUBLIC linalg)
