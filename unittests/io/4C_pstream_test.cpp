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
    CORE::IO::Pstream ps;
    EXPECT_THROW(ps.flush(), CORE::Exception);
    EXPECT_THROW((ps << "blub"), CORE::Exception);
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, DoubleInitializeThrows)
  {
    CORE::IO::Pstream ps;
    ps.setup(true, false, true, CORE::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 4, "");
    EXPECT_THROW(ps.setup(false, false, false, CORE::IO::standard,
                     Teuchos::rcp(new Epetra_SerialComm), 0, 2, ""),
        CORE::Exception);
  }

  TEST(PstreamTest, NonexistantProc)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    EXPECT_THROW(ps.setup(false, false, false, CORE::IO::standard,
                     Teuchos::rcp(new Epetra_SerialComm), 4, 2, ""),
        CORE::Exception);
  }

  TEST(PstreamTest, InitializedUse)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(true, false, false, CORE::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps << "blub");
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, OutputLevel)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(true, false, false, CORE::IO::minimal, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_EQ(ps.requested_output_level(), CORE::IO::minimal);
    CORE::IO::Level &lvl = ps(CORE::IO::debug);
    EXPECT_NO_THROW(lvl << 4);
  }

  TEST(PstreamTest, InputTypes)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(false, false, true, CORE::IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << 4UL << -5LL << 1337.0 << 42.0f << "blub" << std::string("blah") << "\n");
    EXPECT_NO_THROW(ps.flush());
    EXPECT_NO_THROW(ps.close());
  }

  TEST(PstreamTest, ExternalOperators)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(false, false, true, CORE::IO::debug, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    EXPECT_NO_THROW(ps << "blub" << CORE::IO::flush);
    EXPECT_NO_THROW(ps << "blah" << CORE::IO::endl);
  }

  TEST(PstreamTest, Level)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(true, false, true, CORE::IO::undef, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    CORE::IO::Level &lvl = ps(CORE::IO::debug);
    EXPECT_NO_THROW(lvl.stream(1.2));
    EXPECT_NO_THROW(lvl << 4);
    EXPECT_NO_THROW(lvl.SetLevel(CORE::IO::minimal) << 5);
  }

  TEST(PstreamTest, LevelExternalOperators)
  {
    using namespace FourC;
    CORE::IO::Pstream ps;
    ps.setup(true, false, true, CORE::IO::standard, Teuchos::rcp(new Epetra_SerialComm), 0, 0, "");
    CORE::IO::Level &lvl = ps(CORE::IO::debug);
    EXPECT_NO_THROW(lvl << 1.2 << CORE::IO::flush);
    EXPECT_NO_THROW(lvl << 23 << CORE::IO::endl);
  }

}  // namespace