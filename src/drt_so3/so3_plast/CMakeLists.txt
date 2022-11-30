set(OBJS_SOH8_LIB_PLAST
    so3_ssn_plast.cpp
    so3_ssn_plast_evaluate.cpp
    so3_ssn_plast_eletypes.cpp
    so3_ssn_plast_eas.cpp
    so3_ssn_plast_sosh8.cpp
    so3_ssn_plast_sosh18.cpp
    )

baci_add_library(ele_so_plast ${OBJS_SOH8_LIB_PLAST})
target_link_libraries(ele_so_plast PUBLIC ele_so)
