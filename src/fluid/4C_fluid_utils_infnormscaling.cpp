/*-----------------------------------------------------------*/
/*! \file

\brief infnorm-scaling utility class for preconditioning of fluid problems


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_utils_infnormscaling.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparsematrix.hpp"

#include <Epetra_Comm.h>
#include <stdio.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::UTILS::FluidInfNormScaling::FluidInfNormScaling(Core::LinAlg::MapExtractor& mapextractor)
    : myrank_(mapextractor.Map(0)->Comm().MyPID()),
      velpressplitter_(mapextractor),
      leftscale_momentum_(true),
      leftscale_continuity_(false)
{
  return;
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::FluidInfNormScaling::scale_system(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Epetra_Vector& b)
{
  if (myrank_ == 0) std::cout << "Performing scaling of linear system" << std::endl;

  // The matrices are modified here. Do we have to revert the change later on?
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> matrcp =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

  if (matrcp != Teuchos::null)  // yes, we have a block sparse matrix
  {
    Core::LinAlg::BlockSparseMatrixBase& mat = *matrcp;

    Teuchos::RCP<Epetra_CrsMatrix> A00 = mat.matrix(0, 0).epetra_matrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A00->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A00->RowMap(), false));

    if (leftscale_momentum_)
    {
      A00->InvRowSums(*srowsum_);
      if (myrank_ == 0) std::cout << "do left scaling momentum" << std::endl;

      // we want to have the infnorm of the whole(!) row including
      // the off-diagonal block matrix M_(0,1)
      Teuchos::RCP<Epetra_Vector> temp1 =
          Teuchos::rcp(new Epetra_Vector(mat.matrix(0, 0).epetra_matrix()->RowMap(), false));
      srowsum_->Reciprocal(*srowsum_);
      mat.matrix(0, 1).epetra_matrix()->InvRowSums(*temp1);
      temp1->Reciprocal(*temp1);
      srowsum_->Update(1.0, *temp1, 1.0);
      srowsum_->Reciprocal(*srowsum_);
    }
    else
    {
      // no scaling
      srowsum_->PutScalar(1.0);
    }

    scolsum_->PutScalar(1.0);

    if (A00->LeftScale(*srowsum_) or A00->RightScale(*scolsum_) or
        mat.matrix(0, 1).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(1, 0).epetra_matrix()->RightScale(*scolsum_))
      FOUR_C_THROW("fluid scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = velpressplitter_.extract_vector(b, 0);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("fluid scaling failed");

    velpressplitter_.insert_vector(*sx, 0, b);

    // continuity equation
    Teuchos::RCP<Epetra_CrsMatrix> A11 = mat.matrix(1, 1).epetra_matrix();
    prowsum_ = Teuchos::rcp(new Epetra_Vector(A11->RowMap(), false));
    pcolsum_ = Teuchos::rcp(new Epetra_Vector(A11->RowMap(), false));

    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(A11->RowMap(), false));
    if (leftscale_continuity_)
    {
      A11->InvRowSums(*prowsum_);
      if (myrank_ == 0) std::cout << "do left scaling continuity" << std::endl;

      // we want to have the infnorm of the whole(!) row including
      // the off-diagonal block matrix M_(1,0)
      prowsum_->Reciprocal(*prowsum_);
      mat.matrix(1, 0).epetra_matrix()->InvRowSums(*temp);
      temp->Reciprocal(*temp);
      prowsum_->Update(1.0, *temp, 1.0);
      prowsum_->Reciprocal(*prowsum_);
    }
    else
    {
      prowsum_->PutScalar(1.0);
    }

    pcolsum_->PutScalar(1.0);

    if (A11->LeftScale(*prowsum_) or
        //      A->RightScale(*pcolsum_) or
        mat.matrix(1, 0).epetra_matrix()->LeftScale(*prowsum_)  // or
        // mat.Matrix(0,1).EpetraMatrix()->RightScale(*pcolsum_)
    )
      FOUR_C_THROW("fluid scaling failed");

    Teuchos::RCP<Epetra_Vector> px = velpressplitter_.extract_vector(b, 1);

    if (px->Multiply(1.0, *prowsum_, *px, 0.0)) FOUR_C_THROW("fluid scaling failed");

    velpressplitter_.insert_vector(*px, 1, b);


  }  // BlockSparseMatrix

  else  // we have a normal SparseMatrix
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> smat =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(matrix);
    if (smat == Teuchos::null) FOUR_C_THROW("Something went wrong.");

    srowsum_ = Teuchos::rcp(new Epetra_Vector(smat->row_map(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(smat->row_map(), false));
    prowsum_ = Teuchos::null;
    pcolsum_ = Teuchos::null;

    smat->epetra_matrix()->InvRowSums(*srowsum_);
    if (myrank_ == 0) std::cout << "do left scaling of SparseMatrix" << std::endl;

    // leave continuity equation unscaled! -> scaling factors are one
    Teuchos::RCP<Epetra_Vector> px = velpressplitter_.extract_vector(*srowsum_, 1);
    px->PutScalar(1.0);
    velpressplitter_.insert_vector(*px, 1, *srowsum_);

    if (smat->left_scale(*srowsum_)) FOUR_C_THROW("fluid scaling failed");
    if (b.Multiply(1.0, *srowsum_, b, 0.0)) FOUR_C_THROW("fluid scaling failed");

    smat->epetra_matrix()->InvColSums(*scolsum_);
    if (myrank_ == 0) std::cout << "do right scaling pressure" << std::endl;

    // leave velocity columns equation unscaled!
    Teuchos::RCP<Epetra_Vector> ux = velpressplitter_.extract_vector(*scolsum_, 0);
    ux->PutScalar(1.0);
    velpressplitter_.insert_vector(*ux, 0, *scolsum_);

    if (smat->right_scale(*scolsum_)) FOUR_C_THROW("fluid scaling failed");
  }  // SparseMatrix


  // do output
  double srownorm, scolnorm, prownorm = 0.0;
  srowsum_->MeanValue(&srownorm);
  scolsum_->MeanValue(&scolnorm);
  if (prowsum_ != Teuchos::null) prowsum_->MeanValue(&prownorm);

  if (myrank_ == 0)
    std::cout << "MEAN: leftscalemom: " << srownorm << "  rightscale: " << scolnorm
              << "  leftscaleconti: " << prownorm << std::endl;

  srowsum_->MinValue(&srownorm);
  scolsum_->MinValue(&scolnorm);
  if (prowsum_ != Teuchos::null) prowsum_->MinValue(&prownorm);
  if (myrank_ == 0)
    std::cout << "MIN: leftscalemom: " << srownorm << "  rightscale: " << scolnorm
              << "  leftscaleconti: " << prownorm << std::endl;

  srowsum_->MaxValue(&srownorm);
  scolsum_->MaxValue(&scolnorm);
  if (prowsum_ != Teuchos::null) prowsum_->MaxValue(&prownorm);
  if (myrank_ == 0)
    std::cout << "MAX: leftscalemom: " << srownorm << "  rightscale: " << scolnorm
              << "  leftscaleconti: " << prownorm << std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::FluidInfNormScaling::unscale_solution(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Epetra_Vector& x, Epetra_Vector& b)
{
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> matrcp =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

  if (matrcp != Teuchos::null)  // yes, we have a block sparse matrix
  {
    Core::LinAlg::BlockSparseMatrixBase& mat = *matrcp;

    Teuchos::RCP<Epetra_Vector> sy = velpressplitter_.extract_vector(x, 0);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("fluid scaling failed");

    velpressplitter_.insert_vector(*sy, 0, x);

    Teuchos::RCP<Epetra_Vector> sx = velpressplitter_.extract_vector(b, 0);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("fluid scaling failed");

    velpressplitter_.insert_vector(*sx, 0, b);

    Teuchos::RCP<Epetra_CrsMatrix> A00 = mat.matrix(0, 0).epetra_matrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A00->LeftScale(*srowsum_) or A00->RightScale(*scolsum_) or
        mat.matrix(0, 1).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(1, 0).epetra_matrix()->RightScale(*scolsum_))
      FOUR_C_THROW("fluid scaling failed");

    // undo left scaling of continuity equation
    Teuchos::RCP<Epetra_CrsMatrix> A11 = mat.matrix(1, 1).epetra_matrix();
    prowsum_->Reciprocal(*prowsum_);
    if (A11->LeftScale(*prowsum_) or mat.matrix(1, 0).epetra_matrix()->LeftScale(*prowsum_))
      FOUR_C_THROW("fluid scaling failed");
  }
  else
  {
    if (x.Multiply(1.0, *scolsum_, x, 0.0)) FOUR_C_THROW("fluid unscaling failed");

    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);

    // revert matrix and rhs here
    if (myrank_ == 0)
      std::cout << "Only unscaling for solution vector!!! Matrix untouched. " << std::endl;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
