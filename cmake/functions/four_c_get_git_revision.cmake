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

  # get Baci git hash
  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%H
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE BaciGitHash
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(BaciGitHash
        "Unable to determine Baci git hash!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(BaciGitHash
        "${BaciGitHash}"
        PARENT_SCOPE
        )
  endif()

  # get Baci git hash (short)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%h
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE BaciGitHashShort
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(BaciGitHashShort
        "Unable to determine Baci git hash (short)!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(BaciGitHashShort
        "${BaciGitHashShort}"
        PARENT_SCOPE
        )
  endif()

  # get Baci git branch name
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var
    OUTPUT_VARIABLE BaciGitBranch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(NOT ${res_var} EQUAL 0)
    set(BaciGitBranch
        "Unable to determine Baci git branch!"
        PARENT_SCOPE
        )
    message(WARNING "Git command failed! Build will not contain git revision information.")
  else()
    set(BaciGitBranch
        "${BaciGitBranch}"
        PARENT_SCOPE
        )
  endif()

endfunction()
