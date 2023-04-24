set(OBJS_LIB_BEAMCONTACT
    # cmake-format: sortable
    beamcontact_beam3contact.cpp
    beamcontact_beam3contact_manager.cpp
    beamcontact_beam3contact_octtree.cpp
    beamcontact_beam3contactinterface.cpp
    beamcontact_beam3contactnew.cpp
    beamcontact_beam3contactvariables.cpp
    beamcontact_beam3tobeampotential.cpp
    beamcontact_beam3tosolidcontact.cpp
    beamcontact_beam3tosolidmeshtying.cpp
    beamcontact_beam3tospherecontact.cpp
    beamcontact_beam3tospherepotential.cpp
    beamcontact_str_model_evaluator_beaminteraction_old.cpp
    )

baci_add_library(lib_beamcontact ${OBJS_LIB_BEAMCONTACT})

set(LIBS_LIB_BEAMCONTACT
    # cmake-format: sortable
    lib_beaminteraction
    )

baci_add_dependency(lib_beamcontact ${LIBS_LIB_BEAMCONTACT})
