/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_amgnxn_vcycle.hpp"

#include "4C_linear_solver_amgnxn_smoothers.hpp"
#include "4C_utils_exceptions.hpp"

#include <EpetraExt_RowMatrixOut.h>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
CORE::LINEAR_SOLVER::AMGNXN::Vcycle::Vcycle(int NumLevels, int NumSweeps, int FirstLevel)
    : num_levels_(NumLevels),
      num_sweeps_(NumSweeps),
      first_level_(FirstLevel),
      avec_(NumLevels, Teuchos::null),
      pvec_(NumLevels - 1, Teuchos::null),
      rvec_(NumLevels - 1, Teuchos::null),
      svec_pre_(NumLevels, Teuchos::null),
      svec_pos_(NumLevels - 1, Teuchos::null),
      flag_set_up_a_(false),
      flag_set_up_p_(false),
      flag_set_up_r_(false),
      flag_set_up_pre_(false),
      flag_set_up_pos_(false)
{
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::SetOperators(
    std::vector<Teuchos::RCP<BlockedMatrix>> Avec)
{
  if ((int)Avec.size() != num_levels_) FOUR_C_THROW("Error in Setting Avec_: Size dismatch.");
  for (int i = 0; i < num_levels_; i++)
  {
    if (Avec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Avec_: Null pointer.");
    avec_[i] = Avec[i];
  }
  flag_set_up_a_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::SetProjectors(
    std::vector<Teuchos::RCP<BlockedMatrix>> Pvec)
{
  if ((int)Pvec.size() != num_levels_ - 1) FOUR_C_THROW("Error in Setting Pvec_: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (Pvec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Pvec_: Null pointer.");
    pvec_[i] = Pvec[i];
  }
  flag_set_up_p_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::SetRestrictors(
    std::vector<Teuchos::RCP<BlockedMatrix>> Rvec)
{
  if ((int)Rvec.size() != num_levels_ - 1) FOUR_C_THROW("Error in Setting Rvec_: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (Rvec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Rvec_: Null pointer.");
    rvec_[i] = Rvec[i];
  }
  flag_set_up_r_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::SetPreSmoothers(
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPre)
{
  if ((int)SvecPre.size() != num_levels_) FOUR_C_THROW("Error in Setting SvecPre: Size dismatch.");
  for (int i = 0; i < num_levels_; i++)
  {
    if (SvecPre[i] == Teuchos::null) FOUR_C_THROW("Error in Setting SvecPre: Null pointer.");
    svec_pre_[i] = SvecPre[i];
  }
  flag_set_up_pre_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::SetPosSmoothers(
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPos)
{
  if ((int)SvecPos.size() != num_levels_ - 1)
    FOUR_C_THROW("Error in Setting SvecPos: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (SvecPos[i] == Teuchos::null) FOUR_C_THROW("Error in Setting SvecPos: Null pointer.");
    svec_pos_[i] = SvecPos[i];
  }
  flag_set_up_pos_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::DoVcycle(
    const BlockedVector& X, BlockedVector& Y, int level, bool InitialGuessIsZero) const
{
  if (level != num_levels_ - 1)  // Perform one iteration of the V-cycle
  {
    // Apply presmoother
    svec_pre_[level]->Solve(X, Y, InitialGuessIsZero);

    // Compute residual
    BlockedVector DX = X.DeepCopy();
    avec_[level]->Apply(Y, DX);
    DX.Update(1.0, X, -1.0);

    //  Create coarser representation of the residual
    int NV = X.GetVector(0)->NumVectors();
    Teuchos::RCP<BlockedVector> DXcoarse = rvec_[level]->NewRangeBlockedVector(NV, false);
    rvec_[level]->Apply(DX, *DXcoarse);

    // Damp error with coarser levels
    Teuchos::RCP<BlockedVector> DYcoarse = pvec_[level]->NewDomainBlockedVector(NV, false);
    DoVcycle(*DXcoarse, *DYcoarse, level + 1, true);

    // Compute correction
    BlockedVector DY = Y.DeepCopy();
    pvec_[level]->Apply(*DYcoarse, DY);
    Y.Update(1.0, DY, 1.0);

    // Apply post smoother
    svec_pos_[level]->Solve(X, Y, false);
  }
  else  // Apply presmoother
  {
    svec_pre_[level]->Solve(X, Y, InitialGuessIsZero);
  }


  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void CORE::LINEAR_SOLVER::AMGNXN::Vcycle::Solve(
    const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero) const
{
  // Check if everithing is set up
  if (!flag_set_up_a_) FOUR_C_THROW("Operators missing");
  if (!flag_set_up_p_) FOUR_C_THROW("Projectors missing");
  if (!flag_set_up_r_) FOUR_C_THROW("Restrictors missing");
  if (!flag_set_up_pre_) FOUR_C_THROW("Pre-smoothers missing");
  if (!flag_set_up_pos_) FOUR_C_THROW("Post-smoothers missing");

  // Work!
  for (int i = 0; i < num_sweeps_; i++) DoVcycle(X, Y, first_level_, InitialGuessIsZero and i == 0);
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::VcycleSingle(
    int NumLevels, int NumSweeps, int FirstLevel)
    : num_levels_(NumLevels),
      num_sweeps_(NumSweeps),
      first_level_(FirstLevel),
      avec_(NumLevels, Teuchos::null),
      pvec_(NumLevels - 1, Teuchos::null),
      rvec_(NumLevels - 1, Teuchos::null),
      svec_pre_(NumLevels, Teuchos::null),
      svec_pos_(NumLevels - 1, Teuchos::null),
      flag_set_up_a_(false),
      flag_set_up_p_(false),
      flag_set_up_r_(false),
      flag_set_up_pre_(false),
      flag_set_up_pos_(false)
{
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::SetOperators(
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Avec)
{
  if ((int)Avec.size() != num_levels_) FOUR_C_THROW("Error in Setting Avec_: Size dismatch.");
  for (int i = 0; i < num_levels_; i++)
  {
    if (Avec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Avec_: Null pointer.");
    avec_[i] = Avec[i];
  }
  flag_set_up_a_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::SetProjectors(
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Pvec)
{
  if ((int)Pvec.size() != num_levels_ - 1) FOUR_C_THROW("Error in Setting Pvec_: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (Pvec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Pvec_: Null pointer.");
    pvec_[i] = Pvec[i];
  }
  flag_set_up_p_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::SetRestrictors(
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Rvec)
{
  if ((int)Rvec.size() != num_levels_ - 1) FOUR_C_THROW("Error in Setting Rvec_: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (Rvec[i] == Teuchos::null) FOUR_C_THROW("Error in Setting Rvec_: Null pointer.");
    rvec_[i] = Rvec[i];
  }
  flag_set_up_r_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::SetPreSmoothers(
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre)
{
  if ((int)SvecPre.size() != num_levels_) FOUR_C_THROW("Error in Setting SvecPre: Size dismatch.");
  for (int i = 0; i < num_levels_; i++)
  {
    if (SvecPre[i] == Teuchos::null) FOUR_C_THROW("Error in Setting SvecPre: Null pointer.");
    svec_pre_[i] = SvecPre[i];
  }
  flag_set_up_pre_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::SetPosSmoothers(
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos)
{
  if ((int)SvecPos.size() != num_levels_ - 1)
    FOUR_C_THROW("Error in Setting SvecPos: Size dismatch.");
  for (int i = 0; i < num_levels_ - 1; i++)
  {
    if (SvecPos[i] == Teuchos::null) FOUR_C_THROW("Error in Setting SvecPos: Null pointer.");
    svec_pos_[i] = SvecPos[i];
  }
  flag_set_up_pos_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::DoVcycle(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level, bool InitialGuessIsZero) const
{
  if (level != num_levels_ - 1)  // Perform one iteration of the V-cycle
  {
    // Apply presmoother
    svec_pre_[level]->Apply(X, Y, InitialGuessIsZero);

    // Compute residual
    int NV = X.NumVectors();
    Epetra_MultiVector DX(X.Map(), NV, false);
    avec_[level]->Apply(Y, DX);
    DX.Update(1.0, X, -1.0);

    //  Create coarser representation of the residual
    const Epetra_Map& Map = rvec_[level]->RangeMap();
    Epetra_MultiVector DXcoarse(Map, NV, false);
    rvec_[level]->Apply(DX, DXcoarse);

    // Damp error with coarser levels
    const Epetra_Map& Map2 = pvec_[level]->DomainMap();
    Epetra_MultiVector DYcoarse(Map2, NV, false);
    DoVcycle(DXcoarse, DYcoarse, level + 1, true);

    // Compute correction
    Epetra_MultiVector DY(Y.Map(), NV, false);
    pvec_[level]->Apply(DYcoarse, DY);
    Y.Update(1.0, DY, 1.0);

    // Apply post smoother
    svec_pos_[level]->Apply(X, Y, false);
  }
  else  // Apply presmoother
  {
    svec_pre_[level]->Apply(X, Y, InitialGuessIsZero);
  }

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void CORE::LINEAR_SOLVER::AMGNXN::VcycleSingle::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  // Check if everithing is set up
  if (!flag_set_up_a_) FOUR_C_THROW("Operators missing");
  if (!flag_set_up_p_) FOUR_C_THROW("Projectors missing");
  if (!flag_set_up_r_) FOUR_C_THROW("Restrictors missing");
  if (!flag_set_up_pre_) FOUR_C_THROW("Pre-smoothers missing");
  if (!flag_set_up_pos_) FOUR_C_THROW("Post-smoothers missing");

  // Work!
  for (int i = 0; i < num_sweeps_; i++) DoVcycle(X, Y, first_level_, InitialGuessIsZero and i == 0);
  return;
}

FOUR_C_NAMESPACE_CLOSE
