// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UNITTEST_UTILS_SUPPORT_FILES_TEST_HPP
#define FOUR_C_UNITTEST_UTILS_SUPPORT_FILES_TEST_HPP

#include "4C_config.hpp"

#include <filesystem>

namespace TESTING
{
  //! Get the path to the support files for the unit tests.
  inline std::filesystem::path get_support_file_path(std::filesystem::path relative_path)
  {
#ifndef FOUR_C_TEST_SUPPORT_FILE_DIR
#error "FOUR_C_TEST_SUPPORT_FILE_DIR is not defined. Something is wrong in the build system."
#else
    return std::filesystem::path(FOUR_C_TEST_SUPPORT_FILE_DIR) / relative_path;
#endif
  }
}  // namespace TESTING

#endif
