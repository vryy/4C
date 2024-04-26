# get git revision information
function(four_c_get_git_revision_information)

  # check git package
  if(NOT GIT_FOUND)
    find_package(Git QUIET)
  endif()
  if(NOT GIT_FOUND)
    message(WARNING "Git executable not found! Build will not contain git revision information.")
    return()
  endif()

  # enforces reconfigure upon new commit
  if(EXISTS ${PROJECT_SOURCE_DIR}/.git/HEAD)
    configure_file(
      "${PROJECT_SOURCE_DIR}/.git/HEAD" "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/HEAD"
      )
    file(STRINGS ${PROJECT_SOURCE_DIR}/.git/HEAD _head_ref LIMIT_COUNT 1)
    string(REPLACE "ref: " "" _head_ref ${_head_ref})
    if(EXISTS ${PROJECT_SOURCE_DIR}/.git/${_head_ref})
      configure_file(
        "${PROJECT_SOURCE_DIR}/.git/${_head_ref}"
        "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/HEAD_REF"
        )
    endif()
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%H
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE FOUR_C_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(FOUR_C_GIT_HASH
        "Unable to determine 4C git hash!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(FOUR_C_GIT_HASH
        "${FOUR_C_GIT_HASH}"
        PARENT_SCOPE
        )
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%h
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE FOUR_C_GIT_HASH_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(FOUR_C_GIT_HASH_SHORT
        "Unable to determine 4C git hash (short)!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(FOUR_C_GIT_HASH_SHORT
        "${FOUR_C_GIT_HASH_SHORT}"
        PARENT_SCOPE
        )
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE FOUR_C_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(FOUR_C_GIT_BRANCH
        "Unable to determine 4C git branch!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(FOUR_C_GIT_BRANCH
        "${FOUR_C_GIT_BRANCH}"
        PARENT_SCOPE
        )
  endif()

endfunction()
