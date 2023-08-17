# get trilinos version
function(baci_get_trilinos_version)

  # get Trilinos version
  set(TrilinosVersion
      "${Trilinos_VERSION}"
      PARENT_SCOPE
      )

  # get Trilinos git hash
  if(EXISTS "${Trilinos_DIR}/../../../TrilinosRepoVersion.txt")
    file(STRINGS "${Trilinos_DIR}/../../../TrilinosRepoVersion.txt" TrilinosRepoVersionFile)

    list(GET TrilinosRepoVersionFile 1 TrilinosRepoVersionFileLine2)
    separate_arguments(TrilinosRepoVersionFileLine2)
    list(GET TrilinosRepoVersionFileLine2 0 TrilinosSHA)

    set(TrilinosGitHash
        "${TrilinosSHA}"
        PARENT_SCOPE
        )
  else()
    set(TrilinosGitHash
        "Unable to determine Trilinos git hash!"
        PARENT_SCOPE
        )
    message(
      WARNING
        "Trilinos repo version file not found! Build will not contain Trilinos git revision information."
      )
  endif()

endfunction()
