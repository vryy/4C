if(DEFINED ENV{UNITTEST_TIMEOUT})
  set(UNITTEST_TIMEOUT $ENV{UNITTEST_TIMEOUT})
  message("GoogleTest: The unit test timeout is $ENV{UNITTEST_TIMEOUT} s.")
else()
  # default unit test timeout, if not passed as an environment variable
  set(UNITTEST_TIMEOUT 10) # Default unit test timeout

  # Determine scaled global test timeout
  math(EXPR UNITTEST_TIMEOUT "10*${GLOBAL_TEST_TIMEOUT_SCALE}")
  message("The scaled global test timeout is ${UNITTEST_TIMEOUT} s.")

  message(
    "GoogleTest: The unit test timeout is not passed as an environment variable. "
    "It is set to the default value 10 s * ${GLOBAL_TEST_TIMEOUT_SCALE} = ${UNITTEST_TIMEOUT} s."
    )
endif()

include(GoogleTest)
include(cmake/BACIAddGoogleTestExecutable.cmake)
include(cmake/BACILinkGoogleTestNecessaryLibraries.cmake)

# Build the common test functionality
add_subdirectory(common)

# List all test directories here
add_subdirectory(beaminteraction)
add_subdirectory(comm)
add_subdirectory(contact_constitutivelaw)
add_subdirectory(fbi)
add_subdirectory(geometry)
add_subdirectory(geometry_pair)
add_subdirectory(inpar)
add_subdirectory(io)
add_subdirectory(lib)
add_subdirectory(mat)
add_subdirectory(matelast)
add_subdirectory(mixture)
add_subdirectory(mortar)
add_subdirectory(particle_engine)
add_subdirectory(particle_interaction)
add_subdirectory(particle_rigidbody)
add_subdirectory(poromultiphase_scatra)
add_subdirectory(so3)
add_subdirectory(headers)
add_subdirectory(linalg)
add_subdirectory(utils)
