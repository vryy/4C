# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# This function sets the clangcuda mode for ${target} and forwards its configured C++ compiler
# launcher (e.g. ccache) to clangcuda++. Environment variables carry the launcher arguments
# through any intermediate compiler wrapper, such as mpic++, while preserving argument boundaries.
# clangcuda++ applies the launcher to the final Clang compile command.
# Other compiler and rule launchers are left unchanged.
# clangcuda_mode can be either CLANGCUDA_MODE_HOST or CLANGCUDA_MODE_DEVICE
function(set_clangcuda_mode target clangcuda_mode)
  get_target_property(_compiler_launcher ${target} CXX_COMPILER_LAUNCHER)
  if(NOT _compiler_launcher MATCHES "-NOTFOUND$" AND NOT _compiler_launcher STREQUAL "")
    list(LENGTH _compiler_launcher _compiler_launcher_arg_count)

    set(_clangcuda_compiler_launcher
        "${CMAKE_COMMAND}"
        -E
        env
        "CLANGCUDA_COMPILER_LAUNCHER_ARG_COUNT=${_compiler_launcher_arg_count}"
        )

    set(_compiler_launcher_arg_index 0)
    foreach(_compiler_launcher_arg IN LISTS _compiler_launcher)
      list(
        APPEND
        _clangcuda_compiler_launcher
        "CLANGCUDA_COMPILER_LAUNCHER_ARG_${_compiler_launcher_arg_index}=${_compiler_launcher_arg}"
        )
      math(EXPR _compiler_launcher_arg_index "${_compiler_launcher_arg_index} + 1")
    endforeach()

    set_property(TARGET ${target} PROPERTY CXX_COMPILER_LAUNCHER "${_clangcuda_compiler_launcher}")
  endif()

  target_compile_definitions(${target} PRIVATE ${clangcuda_mode})
endfunction()
