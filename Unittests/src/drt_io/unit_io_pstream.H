/*!---------------------------------------------------------------------*

\brief Unittests for the Pstream and Level classes

\level 1

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/

#ifndef UNIT_IO_PSTREAM_H_
#define UNIT_IO_PSTREAM_H_

#include <cxxtest/TestSuite.h>
#include <stdexcept>
#include <Epetra_SerialComm.h>
#include <mpi.h>
#include "src/drt_io/io_pstream.H"

namespace IO
{
  class Pstream_TestSuite;
}

class IO::Pstream_TestSuite : public CxxTest::TestSuite
{
 private:
 public:
  void setUp() {}

  void tearDown() {}

  void test_uninitializedUse()
  {
    Pstream ps;
    TS_ASSERT_THROWS(ps.flush(), std::runtime_error);
    TS_ASSERT_THROWS((ps << "blub"), std::runtime_error);
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void test_doubleInitialize()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 4, "");
    TS_ASSERT_THROWS(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 2, ""),
        std::runtime_error);
  }

  void test_nonexistantProc()
  {
    Pstream ps;
    TS_ASSERT_THROWS(
        ps.setup(false, false, false, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 4, 2, ""),
        std::runtime_error);
  }

  void test_initializedUse()
  {
    Pstream ps;
    ps.setup(true, false, false, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(ps.flush());
    TS_ASSERT_THROWS_NOTHING(ps << "blub");
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void test_outputLevel()
  {
    Pstream ps;
    ps.setup(true, false, false, IO::minimal, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_EQUALS(ps.RequestedOutputLevel(), IO::minimal);
    Level& lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl << 4);
  }

  void test_inputTypes()
  {
    Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(
        ps << 4UL << -5LL << 1337.0 << 42.0f << "blub" << std::string("blah") << "\n");
    TS_ASSERT_THROWS_NOTHING(ps.flush());
    TS_ASSERT_THROWS_NOTHING(ps.close());
  }

  void test_externalOperators()
  {
    Pstream ps;
    ps.setup(false, false, true, IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    TS_ASSERT_THROWS_NOTHING(ps << "blub" << IO::flush);
    TS_ASSERT_THROWS_NOTHING(ps << "blah" << IO::endl);
  }

  void test_Level()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Level& lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl.stream(1.2));
    TS_ASSERT_THROWS_NOTHING(lvl << 4);
    TS_ASSERT_THROWS_NOTHING(lvl.SetLevel(IO::minimal) << 5);
  }

  void test_LevelExternalOperators()
  {
    Pstream ps;
    ps.setup(true, false, true, IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Level& lvl = ps(IO::debug);
    TS_ASSERT_THROWS_NOTHING(lvl << 1.2 << IO::flush);
    TS_ASSERT_THROWS_NOTHING(lvl << 23 << IO::endl);
  }

};  // class Pstream_TestSuite

#endif  // UNIT_IO_PSTREAM_H_
