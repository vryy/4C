/*----------------------------------------------------------------------*/
/*! \file
\brief Unittests for comm_utils implementation
\level 0

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_comm_utils.hpp"

#include "baci_io_pstream.hpp"
#include "baci_lib_utils.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <stdexcept>

BACI_NAMESPACE_OPEN

namespace
{
  Teuchos::RCP<CORE::COMM::Communicators> MockUpCommunicators()
  {
    // mock up for command line to create communicators
    std::vector<std::string> argv{
        "dummyEntryInputFile", "-nptype=separateDatFiles", "-ngroup=2", "-glayout=2,3"};

    return CORE::COMM::CreateComm(argv);
  };

  /**
   * Class to setup parallel vectors which are compared.
   */
  class SetupCompareParallelVectorsTest : public ::testing::Test

  {
   public:
    /**
     * \brief Set up the testing environment.
     */
    SetupCompareParallelVectorsTest()
    {
      communicators_ = MockUpCommunicators();
      IO::cout.setup(false, false, false, IO::standard, communicators_->LocalComm(), 0, 0, "dummy");

      // create arbitrary distributed map within each group
      Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(
          new Epetra_Map(numberOfElementsToDistribute_, 0, *communicators_->LocalComm()));
      epetraVector_ = Teuchos::rcp(new Epetra_Vector(*map, false));

      // fill test Epetra_Vector with entry equals gid
      int numMyEles = map->NumMyElements();
      double* values = new double[numMyEles];
      int* indices = new int[numMyEles];
      for (int lid = 0; lid < numMyEles; ++lid)
      {
        indices[lid] = lid;
        values[lid] = map->GID(lid);
      }
      epetraVector_->ReplaceMyValues(numMyEles, values, indices);
    }

    void TearDown() override { IO::cout.close(); }

   public:
    Teuchos::RCP<CORE::COMM::Communicators> communicators_;
    const int numberOfElementsToDistribute_ = 791;
    Teuchos::RCP<Epetra_Vector> epetraVector_;
  };

  /**
   * Class to setup parallel matrices which are compared.
   */
  class SetupCompareParallelMatricesTest : public ::testing::Test

  {
   public:
    /**
     * \brief Set up the testing environment.
     */
    SetupCompareParallelMatricesTest()
    {
      communicators_ = MockUpCommunicators();
      IO::cout.setup(false, false, false, IO::standard, communicators_->LocalComm(), 0, 0, "dummy");

      // create arbitrary distributed map within each group
      Teuchos::RCP<Epetra_Map> rowmap = Teuchos::rcp(
          new Epetra_Map(numberOfElementsToDistribute_, 0, *communicators_->LocalComm()));
      int approximateNumberOfNonZeroesPerRow = 3;
      epetraCrsMatrix_ = Teuchos::rcp(
          new Epetra_CrsMatrix(::Copy, *rowmap, approximateNumberOfNonZeroesPerRow, false));

      // fill tri-diagonal Epetra_CrsMatrix
      double* values = new double[3];
      int* columnIndices = new int[3];
      int numMyEles = rowmap->NumMyElements();
      for (int lid = 0; lid < numMyEles; ++lid)
      {
        int rowgid = rowmap->GID(lid);
        if (rowgid == 0)  // first global row
        {
          int colIndicesFirstRow[2] = {0, 1};
          double valuesFirstRow[2] = {static_cast<double>(rowgid) + colIndicesFirstRow[0],
              static_cast<double>(rowgid) + colIndicesFirstRow[1]};
          epetraCrsMatrix_->InsertGlobalValues(rowgid, 2, valuesFirstRow, colIndicesFirstRow);
        }
        else if (rowgid == numberOfElementsToDistribute_ - 1)  // last global row
        {
          rowgid = rowmap->GID(numMyEles - 1);
          int colIndicesLastRow[2] = {rowgid - 1, rowgid};
          double valuesLastRow[2] = {static_cast<double>(rowgid) + colIndicesLastRow[0],
              static_cast<double>(rowgid) + colIndicesLastRow[1]};
          epetraCrsMatrix_->InsertGlobalValues(rowgid, 2, valuesLastRow, colIndicesLastRow);
        }
        else  // all rows in between
        {
          columnIndices[0] = rowgid - 1;
          columnIndices[1] = rowgid;
          columnIndices[2] = rowgid + 1;
          values[0] = static_cast<double>(rowgid) + columnIndices[0];
          values[1] = static_cast<double>(rowgid) + columnIndices[1];
          values[2] = static_cast<double>(rowgid) + columnIndices[2];
          epetraCrsMatrix_->InsertGlobalValues(rowgid, 3, values, columnIndices);
        }
      }

      epetraCrsMatrix_->FillComplete(false);
    }

    void TearDown() override { IO::cout.close(); }

   public:
    Teuchos::RCP<CORE::COMM::Communicators> communicators_;
    const int numberOfElementsToDistribute_ = 673;
    Teuchos::RCP<Epetra_CrsMatrix> epetraCrsMatrix_;
  };

  TEST_F(SetupCompareParallelVectorsTest, PositiveTestCompareVectors)
  {
    bool success = CORE::COMM::AreDistributedVectorsIdentical(*communicators_, epetraVector_,  //
        "epetraVector");
    EXPECT_EQ(success, true);
  }

  TEST_F(SetupCompareParallelVectorsTest, NegativeTestCompareVectors)
  {
    // disturb one value on each proc which leads to a failure of the comparison
    const int lastLocalIndex = epetraVector_->MyLength() - 1;
    double disturbedValue = static_cast<double>(lastLocalIndex);
    epetraVector_->ReplaceMyValues(1, &disturbedValue, &lastLocalIndex);

    EXPECT_THROW(
        CORE::COMM::AreDistributedVectorsIdentical(*communicators_, epetraVector_, "epetraVector"),
        CORE::Exception);
  }

  TEST_F(SetupCompareParallelMatricesTest, PositiveTestCompareMatrices)
  {
    bool success = CORE::COMM::AreDistributedSparseMatricesIdentical(
        *communicators_, epetraCrsMatrix_, "epetraCrsMatrix");
    EXPECT_EQ(success, true);
  }

  TEST_F(SetupCompareParallelMatricesTest, NegativeTestCompareMatrices)
  {
    // disturb one value on each proc which leads to a failure of the comparison
    const int myLastLid[1] = {epetraCrsMatrix_->RowMap().NumMyElements() - 1};
    const double value[1] = {static_cast<double>(myLastLid[0])};

    epetraCrsMatrix_->InsertMyValues(myLastLid[0], 1, &value[0], &myLastLid[0]);
    epetraCrsMatrix_->FillComplete(false);

    EXPECT_THROW(CORE::COMM::AreDistributedSparseMatricesIdentical(
                     *communicators_, epetraCrsMatrix_, "epetraCrsMatrix"),
        CORE::Exception);
  }

}  // namespace

BACI_NAMESPACE_CLOSE
