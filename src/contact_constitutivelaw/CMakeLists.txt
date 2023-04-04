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
    contact_constitutivelaw_mirco_contactconstitutivelaw.cpp
    contact_constitutivelaw_power_contactconstitutivelaw.cpp
    )

baci_add_library(lib_contact_constitutivelaw ${OBJS_CONTACT_CONSTITUTIVELAW})

set(LIBS_LIB_CONTACT_CONSTITUTIVELAW
    # cmake-format: sortable
    alg_contact
    )

if(BACI_WITH_MIRCO)
  list(APPEND LIBS_LIB_CONTACT_CONSTITUTIVELAW mirco::mirco_lib)
endif(BACI_WITH_MIRCO)

baci_add_dependency(lib_contact_constitutivelaw ${LIBS_LIB_CONTACT_CONSTITUTIVELAW})
