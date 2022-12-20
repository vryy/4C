set(OBJS_CONTACT_CONSTITUTIVELAW
    # cmake-format: sortable
    brokenrational_contactconstitutivelaw.cpp
    constitutivelaw_definition.cpp
    contact_constitutivelaw_bundle.cpp
    contact_constitutivelaw_interface.cpp
    contactconstitutivelaw.cpp
    contactconstitutivelaw_parameter.cpp
    cubic_contactconstitutivelaw.cpp
    linear_contactconstitutivelaw.cpp
    power_contactconstitutivelaw.cpp
    )

baci_add_library(lib_contact_constitutivelaw ${OBJS_CONTACT_CONSTITUTIVELAW})

set(LIBS_LIB_CONTACT_CONSTITUTIVELAW
    # cmake-format: sortable
    alg_contact
    )

target_link_libraries(lib_contact_constitutivelaw PUBLIC ${LIBS_LIB_CONTACT_CONSTITUTIVELAW})
