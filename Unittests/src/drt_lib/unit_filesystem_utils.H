/*---------------------------------------------------------------------------*/
/*! \file

\brief unittests for utilities concerning file system path operations

\level 3

*/

#ifndef UNIT_FILESYSTEM_UTILS_H_
#define UNIT_FILESYSTEM_UTILS_H_

#include "Unittests/src/common/unit_cxx_test_wrapper.H"
#include <string>
#include <src/drt_lib/filesystem_utils.H>

namespace UTILS
{
  namespace FILESYSTEM
  {
    class FileSystemUtils_TestSuite;
  }
}  // namespace UTILS

class UTILS::FILESYSTEM::FileSystemUtils_TestSuite : public BACICxxTestWrapper
{
 public:
  void Setup() override {}

  void TearDown() override {}

  void TestIsAbsolutePath()
  {
    TS_ASSERT(UTILS::FILESYSTEM::IsAbsolutePath("/abs/path/to/file"))
    TS_ASSERT(!UTILS::FILESYSTEM::IsAbsolutePath("rel/path/to/file"));
  }

  void TestGetDirectoryName()
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
