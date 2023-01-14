set(OBJS_CONTACT_CONSTITUTIVELAW
    # cmake-format: sortable
    contact_constitutivelaw_brokenrational_contactconstitutivelaw.cpp
    contact_constitutivelaw_bundle.cpp
    contact_constitutivelaw_constitutivelaw_definition.cpp
    contact_constitutivelaw_contactconstitutivelaw.cpp
    contact_constitutivelaw_contactconstitutivelaw_parameter.cpp
    contact_constitutivelaw_cubic_contactconstitutivelaw.cpp
    contact_constitutivelaw_interface.cpp
    contact_constitutivelaw_linear_contactconstitutivelaw.cpp
    contact_constitutivelaw_power_contactconstitutivelaw.cpp
    )

baci_add_library(lib_contact_constitutivelaw ${OBJS_CONTACT_CONSTITUTIVELAW})

set(LIBS_LIB_CONTACT_CONSTITUTIVELAW
    # cmake-format: sortable
    alg_contact
    )

target_link_libraries(lib_contact_constitutivelaw PUBLIC ${LIBS_LIB_CONTACT_CONSTITUTIVELAW})
