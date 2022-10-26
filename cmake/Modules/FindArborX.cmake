option(BACI_WITH_ARBORX "Use ArborX for search algorithms" OFF)
if(BACI_WITH_ARBORX)

  message(STATUS "Search algorithms with ArborX: enabled")

  # Unconditionally turn on MPI support inside ArborX
  set(ARBORX_ENABLE_MPI "ON")
  fetchcontent_declare(
    arborx
    GIT_REPOSITORY https://github.com/arborx/ArborX.git
    GIT_TAG v1.2
    )
  fetchcontent_makeavailable(arborx)

  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS ArborX::ArborX)

  add_definitions("-DHAVE_ARBORX")
  set(HAVE_ARBORX ON)
else()
  message(STATUS "Search algorithms with ArborX: disabled")
endif()
