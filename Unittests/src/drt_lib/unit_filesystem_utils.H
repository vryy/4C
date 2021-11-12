/*---------------------------------------------------------------------------*/
/*! \file

\brief unittests for utilities concerning file system path operations

\level 3

*/

#ifndef UNIT_FILESYSTEM_UTILS_H_
#define UNIT_FILESYSTEM_UTILS_H_

#include <cxxtest/TestSuite.h>
#include <string>
#include <src/drt_lib/filesystem_utils.H>

namespace UTILS
{
  namespace FILESYSTEM
  {
    class FileSystemUtils_TestSuite;
  }
}  // namespace UTILS

class UTILS::FILESYSTEM::FileSystemUtils_TestSuite : public CxxTest::TestSuite
{
 public:
  void setUp() override {}

  void tearDown() override {}

  void test_IsAbsolutePath()
  {
    TS_ASSERT(UTILS::FILESYSTEM::IsAbsolutePath("/abs/path/to/file"))
    TS_ASSERT(!UTILS::FILESYSTEM::IsAbsolutePath("rel/path/to/file"))
  }

  void test_GetDirectoryName()
  {
    TS_ASSERT(UTILS::FILESYSTEM::GetDirectoryName("/abs/path/to/file") == "/abs/path/to/")
    TS_ASSERT(UTILS::FILESYSTEM::GetDirectoryName("file") == "./")
    TS_ASSERT(UTILS::FILESYSTEM::GetDirectoryName("rel/path/to/file") == "rel/path/to/")
  }

  void test_Join()
  {
    TS_ASSERT(UTILS::FILESYSTEM::Join("path/to", "file") == "path/to/file")
    TS_ASSERT(UTILS::FILESYSTEM::Join("path/to/", "file") == "path/to/file")
    TS_ASSERT(UTILS::FILESYSTEM::Join("path/to/", "/file") == "/file")
  }
};

#endif
