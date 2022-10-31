if(DOXYGEN_USE_MATHJAX)
  message(STATUS "DOXYGEN: Use mathjax for math equations")

  if(DOXYGEN_USE_LOCAL_MATHJAX)
    set(DOXYGEN_MATHJAX_RELPATH ./mathjax)

    # find local mathjax
    find_file(
      DOXYGEN_MATHJAX_MAIN_FILE MathJax.js
      PATHS ${DOXYGEN_LOCAL_MATHJAX_BASEPATH}
      NO_DEFAULT_PATH
      )

    if(DOXYGEN_MATHJAX_MAIN_FILE)
      message(STATUS "DOXYGEN: MathJax.js found at ${DOXYGEN_LOCAL_MATHJAX_BASEPATH}")
    else()
      message(
        FATAL_ERROR "DOXYGEN: MathJax.js could not be found at ${DOXYGEN_LOCAL_MATHJAX_BASEPATH}"
        )
    endif()
  else()
    set(DOXYGEN_MATHJAX_RELPATH https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/)
    message(STATUS "DOXYGEN: Use Mathjax from CDN")
  endif()
else()
  message(STATUS "DOXYGEN: Use latex for math equations")
endif()
