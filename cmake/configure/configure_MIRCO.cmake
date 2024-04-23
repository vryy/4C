option(FOUR_C_MIRCO_FIND_INSTALLED "Use installed MIRCO instead of fetching sources" OFF)
if(FOUR_C_MIRCO_FIND_INSTALLED)

  message(STATUS "FOUR_C_MIRCO_FIND_INSTALLED is enabled")

  # MIRCO provides a package configuration file if installed.
  find_package(mirco_lib HINTS ${FOUR_C_MIRCO_ROOT})

  if(NOT mirco_lib_FOUND)
    message(
      FATAL_ERROR
        "mirco_lib could not be found. Please ensure that the FOUR_C_MIRCO_ROOT path is correctly defined in the config file. Also, please use 'make install' and not just 'make' to install MIRCO."
      )
  endif()

else() # Fetch MIRCO from GIT repository
  # Turn off googletest and Trilinos in MIRCO so that they don't interfere with BACI
  set(GTEST_IN_MIRCO "OFF")
  set(TRILINOS_IN_MIRCO "OFF")

  fetchcontent_declare(
    mirco
    GIT_REPOSITORY https://github.com/imcs-compsim/MIRCO.git
    GIT_TAG 8a8ae9c703a762459995d56d36292b77ff6c1985
    )
  fetchcontent_makeavailable(mirco)
endif()
