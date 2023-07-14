include(FindPackageHandleStandardArgs)

find_package(Python COMPONENTS Interpreter)

execute_process(COMMAND ${Python_EXECUTABLE} -m venv "${BACI_READTHEDOCS_VENV_DIR}")
set(ENV{VIRTUAL_ENV} "${BACI_READTHEDOCS_VENV_DIR}")
# We are likely to find Sphinx near the Python interpreter
set(Python_FIND_VIRTUALENV ONLY)
unset(Python_EXECUTABLE)
## Launch a new search
find_package(Python COMPONENTS Interpreter Development)

execute_process(
  COMMAND
    ${Python_EXECUTABLE} -m pip install -r ${BACI_READTHEDOCS_CONFIG_DIR}/requirements.txt
    --upgrade pip
  )

find_program(
  SPHINX_EXECUTABLE
  NAMES "${BACI_READTHEDOCS_VENV_DIR}/bin/sphinx-multibuild"
  DOC "Sphinx documentation generator"
  )
message("Found sphinx-multibuild as ${SPHINX_EXECUTABLE}")

mark_as_advanced(SPHINX_EXECUTABLE)

find_package_handle_standard_args(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)
