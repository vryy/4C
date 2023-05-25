baci_auto_define_module()

set(_contact_constitutivelaw_dependencies
    # cmake-format: sortable
    contact
    inpar
    lib
    mat
    )

if(BACI_WITH_MIRCO)
  list(APPEND _contact_constitutivelaw_dependencies mirco::mirco_lib)
endif(BACI_WITH_MIRCO)

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_contact_constitutivelaw_dependencies})
