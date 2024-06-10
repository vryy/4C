#
# Detect which linkers are available and use the fastest one.
#

# The following options may be disabled by users. In this case, we do not search for a linker and use the user flags.
option(FOUR_C_DETECT_LINKER "Detect a fast linker" ON)

if(FOUR_C_DETECT_LINKER)
  # The order of the linkers here is the order in which we prefer to use the linkers.
  # As soon as a linker works, skip the rest.
  set(_linkers "mold" "lld" "gold" "bfd")
  set(_no_linker_found True)
  foreach(_linker_name ${_linkers})
    # Check that the linker exists
    find_program(FOUR_C_HAVE_LINKER_PROGRAM_${_linker_name} "ld.${_linker_name}")

    if(FOUR_C_HAVE_LINKER_PROGRAM_${_linker_name})
      # Check if the linker works out of the box.
      four_c_check_compiles(
        FOUR_C_LINKER_FUNCTIONAL_${_linker_name}
        LINK_OPTIONS
        "-fuse-ld=${_linker_name}"
        APPEND_ON_SUCCESS
        )

      if(FOUR_C_LINKER_FUNCTIONAL_${_linker_name})
        message(STATUS "Linker ${_linker_name} is functional and will be used.")
        set(_no_linker_found False)
        break()
      else()
        # This is a special case which happens when using the mpic++ compiler wrapper on
        # Ubuntu 20.04. The faster linkers can fail due to a mistake in OpenMPI.
        # Since we know how to fix linking, we try if the linker works when adding the missing -lopen-pal flag.
        four_c_check_compiles(
          FOUR_C_LINKER_FUNCTIONAL_WITH_OPEN_PAL_${_linker_name}
          LINK_LIBRARIES
          "open-pal"
          LINK_OPTIONS
          "-fuse-ld=${_linker_name}"
          APPEND_ON_SUCCESS
          )

        if(FOUR_C_LINKER_FUNCTIONAL_WITH_OPEN_PAL_${_linker_name})
          message(
            STATUS
              "Linker ${_linker_name} is functional (when additionally adding `-lopen-pal`) and will be used."
            )
          set(_no_linker_found False)
          break()
        else()
          # Something else must be wrong but we have no more heuristics to fix the problem.
          message(STATUS "Linker ${_linker_name} not functional. Trying another linker.")
        endif()
      endif()
    else()
      message(STATUS "Linker ${_linker_name} not installed. Trying another linker.")
    endif()
  endforeach()

  # This really should not happen since bfd is always available but it does not hurt to check.
  if(_no_linker_found)
    message(
      FATAL_ERROR
        "Failed to find any working linker. Please check your compiler and any manually added flags."
      )
  endif()
endif()
