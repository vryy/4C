set(OBJS_SOH8_LIB_PLAST
    # cmake-format: sortable
    so3_ssn_plast.cpp
    so3_ssn_plast_eas.cpp
    so3_ssn_plast_eletypes.cpp
    so3_ssn_plast_evaluate.cpp
    so3_ssn_plast_sosh18.cpp
    so3_ssn_plast_sosh8.cpp
    )

baci_add_library(ele_so_plast ${OBJS_SOH8_LIB_PLAST})

set(LIBS_ELE_SO_PLAST
    # cmake-format: sortable
    drt_lib
    ele_so
    )

target_link_libraries(ele_so_plast PUBLIC ${LIBS_ELE_SO_PLAST})
