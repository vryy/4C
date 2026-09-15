# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(C4_LIBRARY_TYPE
    "STATIC"
    CACHE STRING "" FORCE
    )

message(STATUS "Fetch content for ryml")
set(C4_LIBRARY_TYPE
    "STATIC"
    CACHE STRING "" FORCE
    )
fetchcontent_declare(
  ryml
  GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
  GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311 # version 0.9.0
  )
set(RYML_INSTALL
    ON
    CACHE BOOL "Turn on ryml install" FORCE
    )
if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(MSVC_VERSION)
    add_compile_definitions(C4_MSVC=1 _MSC_VER=${MSVC_VERSION})
  endif()
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "MSVC") # for Clang-cl
    add_compile_options(/clang:-Wno-nan-infinity-disabled)
    add_compile_options(/clang:-Wno-c++20-extensions)
  endif()
endif()
fetchcontent_makeavailable(ryml)
set_target_properties(ryml c4core PROPERTIES POSITION_INDEPENDENT_CODE ON)
set(FOUR_C_RYML_ROOT "${CMAKE_INSTALL_PREFIX}")

four_c_add_external_dependency(four_c_all_enabled_external_dependencies ryml::ryml)

four_c_remember_variable_for_install(FOUR_C_RYML_ROOT)
