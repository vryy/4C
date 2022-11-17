set(OBJS_DRT_CONTACT_CONSTITUTIVELAW
    contactconstitutivelaw.cpp
    contactconstitutivelaw_parameter.cpp
    brokenrational_contactconstitutivelaw.cpp
    cubic_contactconstitutivelaw.cpp
    linear_contactconstitutivelaw.cpp
    power_contactconstitutivelaw.cpp
    contact_constitutivelaw_bundle.cpp
    contact_constitutivelaw_interface.cpp
    )

baci_add_library(drt_lib_contact_constitutivelaw ${OBJS_DRT_CONTACT_CONSTITUTIVELAW})

target_link_libraries(drt_lib_contact_constitutivelaw PUBLIC alg_contact)
