if(DEFINED ENV{UNITTEST_TIMEOUT})
  set(UNITTEST_TIMEOUT $ENV{UNITTEST_TIMEOUT})
  message("GoogleTest: The unit test timeout is $ENV{UNITTEST_TIMEOUT} s.")
else()
  # default unit test timeout, if not passed as an environment variable
  set(UNITTEST_TIMEOUT 10) # Default unit test timeout

  message(
    "GoogleTest: The unit test timeout is not passed as an environment variable. "
    "It is set to the default value ${UNITTEST_TIMEOUT} s."
    )
endif()

include(GoogleTest)
include(cmake/BACIAddGoogleTestExecutable.cmake)
include(cmake/BACILinkGoogleTestNecessaryLibraries.cmake)

# Build the common test functionality
add_subdirectory(common)

# List all test directories here
add_subdirectory(drt_beaminteraction)
add_subdirectory(drt_contact_constitutivelaw)
add_subdirectory(drt_fbi)
add_subdirectory(drt_geometry_pair)
add_subdirectory(drt_inpar)
add_subdirectory(drt_io)
add_subdirectory(drt_lib)
add_subdirectory(drt_mat)
add_subdirectory(drt_matelast)
add_subdirectory(drt_mortar)
add_subdirectory(drt_particle_engine)
add_subdirectory(drt_particle_interaction)
add_subdirectory(drt_particle_rigidbody)
add_subdirectory(drt_poromultiphase_scatra)
add_subdirectory(drt_so3)
add_subdirectory(headers)
add_subdirectory(linalg)
