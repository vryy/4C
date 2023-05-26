##
# Process global options so they are added to BACI_GLOBAL_COMPILE_DEFINITIONS
#
baci_process_global_option(
  TRILINOS_DEVELOP
  "Select Trilinos installation based on current develop branch (highly experimental!)"
  OFF
  )
baci_process_global_option(DSERROR_DUMP "dserror creates a core file" OFF)
baci_process_global_option(TRAP_FE "Crash BACI if a nan or inf occurs" ON)
