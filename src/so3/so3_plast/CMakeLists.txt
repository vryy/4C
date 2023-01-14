set(OBJS_SOH8_LIB_PLAST
    # cmake-format: sortable
    so3_plast_ssn.cpp
    so3_plast_ssn_eas.cpp
    so3_plast_ssn_eletypes.cpp
    so3_plast_ssn_evaluate.cpp
    so3_plast_ssn_sosh18.cpp
    so3_plast_ssn_sosh8.cpp
    )

baci_add_library(ele_so_plast ${OBJS_SOH8_LIB_PLAST})

set(LIBS_ELE_SO_PLAST
    # cmake-format: sortable
    ele_so
    lib
    )

target_link_libraries(ele_so_plast PUBLIC ${LIBS_ELE_SO_PLAST})
