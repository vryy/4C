/*!----------------------------------------------------------------------
\file teko_baciepetraoperatorwrapper.cpp

\brief wrapper class around the epetra operator in baci for use with teko
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifdef HAVE_TEKO

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"

// Thyra-Epetra adapter includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_EpetraHelpers.hpp"


#include "teko_baciepetraoperatorwrapper.H"


LINALG::SOLVER::TEKO::Teko_BACIMappingStrategy::Teko_BACIMappingStrategy(
    const LINALG::MultiMapExtractor& rangemaps, const LINALG::MultiMapExtractor& domainmaps)
    : rangemaps_(rangemaps), domainmaps_(domainmaps)
{
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      X       - source Epetra_MultiVector
//      thyra_X - destination Thyra::MultiVectorBase
//
void LINALG::SOLVER::TEKO::Teko_BACIMappingStrategy::copyEpetraIntoThyra(
    const Epetra_MultiVector& X, const Teuchos::Ptr<Thyra::MultiVectorBase<double>>& thyra_X) const
{
  Teuchos::Ptr<Thyra::ProductMultiVectorBase<double>> prod_X =
      Teuchos::ptr_dynamic_cast<Thyra::ProductMultiVectorBase<double>>(thyra_X);

  for (int i = 0; i < rangemaps_.NumMaps(); i++)
  {
    Teuchos::RCP<Thyra::DefaultSpmdMultiVector<double>> vec =
        Teuchos::rcp_dynamic_cast<Thyra::DefaultSpmdMultiVector<double>>(
            prod_X->getNonconstMultiVectorBlock(i));
    Teuchos::RCP<Epetra_MultiVector> ep_subX = rangemaps_.ExtractVector(X, i);
    Teko::Epetra::fillDefaultSpmdMultiVector(vec, ep_subX);
  }
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      thyra_Y - source Thyra::MultiVectorBase
//      Y       - destination Epetra_MultiVector
//
void LINALG::SOLVER::TEKO::Teko_BACIMappingStrategy::copyThyraIntoEpetra(
    const Teuchos::RCP<const Thyra::MultiVectorBase<double>>& thyra_Y, Epetra_MultiVector& Y) const
{
  Teuchos::RCP<const Thyra::DefaultProductMultiVector<double>> prod_Y =
      Teuchos::rcp_dynamic_cast<const Thyra::DefaultProductMultiVector<double>>(thyra_Y);

  for (int i = 0; i < rangemaps_.NumMaps(); i++)
  {
    Teuchos::RCP<const Epetra_MultiVector> ep_subY =
        Thyra::get_Epetra_MultiVector(*(rangemaps_.Map(i)), prod_Y->getMultiVectorBlock(i));
    rangemaps_.InsertVector(*ep_subY, i, Y);
  }
}

LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper::Teko_BACIEpetraOperatorWrapper(
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase>& mat)
    : Teko::Epetra::EpetraOperatorWrapper(), mat_(mat)
{
  mappingStrategy_ = Teuchos::rcp(new LINALG::SOLVER::TEKO::Teko_BACIMappingStrategy(
      mat->RangeExtractor(), mat->DomainExtractor()));
  SetMapStrategy(mappingStrategy_);
  BuildBlockedOperator();
  comm_ = getEpetraComm(*thyraOp_);
}

void LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper::BuildBlockedOperator()
{
  /*Teuchos::RCP<const Thyra::LinearOpBase<double> > thA11 =
  Thyra::epetraLinearOp(mat_->Matrix(0,0).EpetraMatrix(),"A11"); Teuchos::RCP<const
  Thyra::LinearOpBase<double> > thA12 =
  Thyra::epetraLinearOp(mat_->Matrix(0,1).EpetraMatrix(),"A12"); Teuchos::RCP<const
  Thyra::LinearOpBase<double> > thA21 =
  Thyra::epetraLinearOp(mat_->Matrix(1,0).EpetraMatrix(),"A21"); Teuchos::RCP<const
  Thyra::LinearOpBase<double> > thA22 =
  Thyra::epetraLinearOp(mat_->Matrix(1,1).EpetraMatrix(),"A22");*/

  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blockA =
      Thyra::defaultBlockedLinearOp<double>();

  int nRows = mat_->Rows();
  int nCols = mat_->Cols();

  blockA->beginBlockFill(nRows, nCols);

  for (int i = 0; i < nRows; i++)
  {
    for (int j = 0; j < nCols; j++)
    {
      std::stringstream label;
      label << "A" << i << j;
      blockA->setBlock(i, j, Thyra::epetraLinearOp(mat_->Matrix(i, j).EpetraMatrix(), label.str()));
    }
  }

  blockA->endBlockFill();
  blockA->setObjectLabel("blockA");

  // Teuchos::RCP<const Thyra::LinearOpBase<double> > blockA =
  // Thyra::block2x2(thA11,thA12,thA21,thA22,"A");
  Teko::BlockedLinearOp A = Teko::toBlockedLinearOp(blockA);
  thyraOp_ = A;
}

Teuchos::RCP<const Thyra::EpetraLinearOp>
LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper::GetThyraBlock(int r, int c)
{
  const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double>> blkOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double>>(getThyraOp());
  return Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(blkOp->getBlock(r, c));
}

#endif /* HAVE_TEKO */
