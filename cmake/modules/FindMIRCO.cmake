# Find the MIRCO includes and libraries
option(BACI_WITH_MIRCO "Use MIRCO as a contact constitutive law" OFF)
if(BACI_WITH_MIRCO)

  message(STATUS "MIRCO is: enabled")
  add_definitions("-DBACI_WITH_MIRCO")

  option(MIRCO_DEVELOP "Use MIRCO experimental" OFF)
  if(MIRCO_DEVELOP)

    message(STATUS "MIRCO_DEVELOP is: enabled")

    find_package(mirco_lib HINTS ${MIRCO_DEVELOP_INSTALL_DIR}/lib/cmake/mirco)

    if(NOT mirco_lib_FOUND)
      message(
        FATAL_ERROR
          "mirco_lib could not be found. Please ensure that the MIRCO_DEVELOP_INSTALL_DIR path is correctly defined in the config file. Also, please use 'make install' and not just 'make' to install MIRCO."
        )
    endif()

  else() # Fetch MIRCO from GIT repository

    # Turn off googletest and Trilinos in MIRCO so that they don't interfere with BACI
    set(GTEST_IN_MIRCO "OFF")
    set(TRILINOS_IN_MIRCO "OFF")

    fetchcontent_declare(
      mirco
      GIT_REPOSITORY https://github.com/imcs-compsim/MIRCO.git
      GIT_TAG 03a8ab0129d9c75748036418230d1119441aae79 #v0.1.0
      )
    fetchcontent_makeavailable(mirco)

  endif()

else()
  message(STATUS "MIRCO is: disabled")
endif(BACI_WITH_MIRCO)
