set(OBJS_DRT_LIB_BEAMCONTACT
    # cmake-format: sortable
    beam3contact.cpp
    beam3contact_manager.cpp
    beam3contact_octtree.cpp
    beam3contactinterface.cpp
    beam3contactnew.cpp
    beam3contactvariables.cpp
    beam3tobeampotential.cpp
    beam3tosolidcontact.cpp
    beam3tosolidmeshtying.cpp
    beam3tospherecontact.cpp
    beam3tospherepotential.cpp
    str_model_evaluator_beaminteraction_old.cpp
    )

baci_add_library(drt_lib_beamcontact ${OBJS_DRT_LIB_BEAMCONTACT})

set(LIBS_DRT_LIB_BEAMCONTACT
    # cmake-format: sortable
    drt_lib_beaminteraction
    )

target_link_libraries(drt_lib_beamcontact PUBLIC ${LIBS_DRT_LIB_BEAMCONTACT})
