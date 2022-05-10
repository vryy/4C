/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the Pstream and Level classes

\level 1

*-----------------------------------------------------------------------*/

#ifndef UNIT_IO_PSTREAM_H_
#define UNIT_IO_PSTREAM_H_

#include "src/common/unit_cxx_test_wrapper.H"
#include <stdexcept>
#include <Epetra_SerialComm.h>
#include <mpi.h>
#include "src/drt_io/io_pstream.H"

namespace IO
{
  class Pstream_TestSuite;
}

class IO::Pstream_TestSuite : public BACICxxTestWrapper
{
 private:
 public:
  void Setup() override {}

  void TearDown() override {}

  void Test_uninitializedUse()
  {
    Pstream ps;
    TS_ASSERT_THROWS(ps.flush(), const std::runtime_error &);
    TS_ASSERT_THROWS((ps << "blub"), const std::runtime_error &);
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void TestDoubleInitialize()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 4, "");
    TS_ASSERT_THROWS(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 2, ""),
        const std::runtime_error &);
  }

  void TestNonexistantProc()
  {
    Pstream ps;
    TS_ASSERT_THROWS(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 4, 2, ""),
        const std::runtime_error &);
  }

  void TestInitializedUse()
  {
    Pstream ps;
    ps.setup(true, false, false, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(ps.flush());
    TS_ASSERT_THROWS_NOTHING(ps << "blub");
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void TestOutputLevel()
  {
    Pstream ps;
    ps.setup(true, false, false, IO::minimal, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_EQUALS(ps.RequestedOutputLevel(), IO::minimal);
    Level &lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl << 4);
  }

  void TestInputTypes()
  {
    Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(
        ps << 4UL << -5LL << 1337.0 << 42.0f << "blub" << std::string("blah") << "\n");
    TS_ASSERT_THROWS_NOTHING(ps.flush());
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void TestExternalOperators()
  {
    Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(ps << "blub" << IO::flush);
    TS_ASSERT_THROWS_NOTHING(ps << "blah" << IO::endl);
  }

  void TestLevel()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Level &lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl.stream(1.2));
    TS_ASSERT_THROWS_NOTHING(lvl << 4);
    TS_ASSERT_THROWS_NOTHING(lvl.SetLevel(IO::minimal) << 5);
  }

  void TestLevelExternalOperators()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Level &lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl << 1.2 << IO::flush);
    TS_ASSERT_THROWS_NOTHING(lvl << 23 << IO::endl);
  }

};  // class Pstream_TestSuite

#endif  // UNIT_IO_PSTREAM_H_
