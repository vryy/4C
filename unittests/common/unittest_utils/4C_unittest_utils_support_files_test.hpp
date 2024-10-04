/*----------------------------------------------------------------------*/
/*! \file
\brief Get the path to the support files for the unit tests
\level 1
*----------------------------------------------------------------------*/
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
