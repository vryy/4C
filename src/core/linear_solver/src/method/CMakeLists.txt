four_c_auto_define_module()

# Clashing definitions of `View` inside Trilinos
four_c_unity_build_compile_separately(
  ${AUTO_DEFINED_MODULE_NAME} 4C_linear_solver_method_direct.cpp
  )
