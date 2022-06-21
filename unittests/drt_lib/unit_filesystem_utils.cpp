/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for utilities concerning file system path operations

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <string>
#include "src/drt_lib/filesystem_utils.H"

namespace
{
  TEST(FilesystemUtilsTest, IsAbsolutePath)
  {
    EXPECT_TRUE(UTILS::FILESYSTEM::IsAbsolutePath("/abs/path/to/file"));
    EXPECT_TRUE(!UTILS::FILESYSTEM::IsAbsolutePath("rel/path/to/file"));
  }

  TEST(FilesystemUtilsTest, GetDirectoryName)
  {
    EXPECT_TRUE(UTILS::FILESYSTEM::GetDirectoryName("/abs/path/to/file") == "/abs/path/to/");
    EXPECT_TRUE(UTILS::FILESYSTEM::GetDirectoryName("file") == "./");
    EXPECT_TRUE(UTILS::FILESYSTEM::GetDirectoryName("rel/path/to/file") == "rel/path/to/");
  };

  TEST(FilesystemUtilsTest, _Join)
  {
    EXPECT_TRUE(UTILS::FILESYSTEM::Join("path/to", "file") == "path/to/file");
    EXPECT_TRUE(UTILS::FILESYSTEM::Join("path/to/", "file") == "path/to/file");
    EXPECT_TRUE(UTILS::FILESYSTEM::Join("path/to/", "/file") == "/file");
  }
}  // namespace