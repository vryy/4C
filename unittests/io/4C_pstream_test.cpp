/*----------------------------------------------------------------------*/
/*! \file
\brief Unittests for the Pstream and Level classes
\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_io_pstream.hpp"

#include <Epetra_SerialComm.h>

#include <stdexcept>

namespace
{
  using namespace FourC;

  TEST(PstreamTest, UninitializedUseThrows)
  {
    Core::IO::Pstream ps;
    EXPECT_THROW(ps.flush(), Core::Exception);
    EXPECT_THROW((ps << "blub"), Core::Exception);
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, DoubleInitializeThrows)
  {
    Core::IO::Pstream ps;
    ps.setup(true, false, true, Core::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 4, "");
    EXPECT_THROW(ps.setup(false, false, false, Core::IO::standard,
                     Teuchos::rcp(new Epetra_SerialComm), 0, 2, ""),
        Core::Exception);
  }

  TEST(PstreamTest, NonexistantProc)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    EXPECT_THROW(ps.setup(false, false, false, Core::IO::standard,
                     Teuchos::rcp(new Epetra_SerialComm), 4, 2, ""),
        Core::Exception);
  }

  TEST(PstreamTest, InitializedUse)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(true, false, false, Core::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps << "blub");
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, OutputLevel)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(true, false, false, Core::IO::minimal, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_EQ(ps.requested_output_level(), Core::IO::minimal);
    Core::IO::Level &lvl = ps(Core::IO::debug);
    EXPECT_NO_THROW(lvl << 4);
  }

  TEST(PstreamTest, InputTypes)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(false, false, true, Core::IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << 4UL << -5LL << 1337.0 << 42.0f << "blub" << std::string("blah") << "\n");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, ExternalOperators)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(false, false, true, Core::IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << "blub" << Core::IO::flush);
    EXPECT_NO_THROW(ps << "blah" << Core::IO::endl);
  }

  TEST(PstreamTest, Level)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(true, false, true, Core::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Core::IO::Level &lvl = ps(Core::IO::debug);
    EXPECT_NO_THROW(lvl.stream(1.2));
    EXPECT_NO_THROW(lvl << 4);
    EXPECT_NO_THROW(lvl.SetLevel(Core::IO::minimal) << 5);
  }

  TEST(PstreamTest, LevelExternalOperators)
  {
    using namespace FourC;
    Core::IO::Pstream ps;
    ps.setup(true, false, true, Core::IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    Core::IO::Level &lvl = ps(Core::IO::debug);
    EXPECT_NO_THROW(lvl << 1.2 << Core::IO::flush);
    EXPECT_NO_THROW(lvl << 23 << Core::IO::endl);
  }

}  // namespace