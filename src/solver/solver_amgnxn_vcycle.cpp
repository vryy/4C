/*!----------------------------------------------------------------------
\file solver_amgnxn_vcycle.cpp

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/


#ifdef HAVE_MueLu

#include <iostream>

#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_EpetraOperator.hpp>
#include "EpetraExt_RowMatrixOut.h"
#include "../drt_lib/drt_dserror.H"
#include "solver_amgnxn_vcycle.H"
#include "solver_amgnxn_smoothers.H"



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMGNXN::Vcycle::Vcycle(int NumLevels, int NumSweeps, int FirstLevel)
    : NumLevels_(NumLevels),
      NumSweeps_(NumSweeps),
      FirstLevel_(FirstLevel),
      Avec_(NumLevels, Teuchos::null),
      Pvec_(NumLevels - 1, Teuchos::null),
      Rvec_(NumLevels - 1, Teuchos::null),
      SvecPre_(NumLevels, Teuchos::null),
      SvecPos_(NumLevels - 1, Teuchos::null),
      flag_set_up_A_(false),
      flag_set_up_P_(false),
      flag_set_up_R_(false),
      flag_set_up_Pre_(false),
      flag_set_up_Pos_(false)
{
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Vcycle::SetOperators(std::vector<Teuchos::RCP<BlockedMatrix>> Avec)
{
  if ((int)Avec.size() != NumLevels_) dserror("Error in Setting Avec_: Size dismatch.");
  for (int i = 0; i < NumLevels_; i++)
  {
    if (Avec[i] == Teuchos::null) dserror("Error in Setting Avec_: Null pointer.");
    Avec_[i] = Avec[i];
  }
  flag_set_up_A_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Vcycle::SetProjectors(std::vector<Teuchos::RCP<BlockedMatrix>> Pvec)
{
  if ((int)Pvec.size() != NumLevels_ - 1) dserror("Error in Setting Pvec_: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (Pvec[i] == Teuchos::null) dserror("Error in Setting Pvec_: Null pointer.");
    Pvec_[i] = Pvec[i];
  }
  flag_set_up_P_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Vcycle::SetRestrictors(std::vector<Teuchos::RCP<BlockedMatrix>> Rvec)
{
  if ((int)Rvec.size() != NumLevels_ - 1) dserror("Error in Setting Rvec_: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (Rvec[i] == Teuchos::null) dserror("Error in Setting Rvec_: Null pointer.");
    Rvec_[i] = Rvec[i];
  }
  flag_set_up_R_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Vcycle::SetPreSmoothers(
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPre)
{
  if ((int)SvecPre.size() != NumLevels_) dserror("Error in Setting SvecPre: Size dismatch.");
  for (int i = 0; i < NumLevels_; i++)
  {
    if (SvecPre[i] == Teuchos::null) dserror("Error in Setting SvecPre: Null pointer.");
    SvecPre_[i] = SvecPre[i];
  }
  flag_set_up_Pre_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Vcycle::SetPosSmoothers(
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPos)
{
  if ((int)SvecPos.size() != NumLevels_ - 1) dserror("Error in Setting SvecPos: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (SvecPos[i] == Teuchos::null) dserror("Error in Setting SvecPos: Null pointer.");
    SvecPos_[i] = SvecPos[i];
  }
  flag_set_up_Pos_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMGNXN::Vcycle::DoVcycle(
    const BlockedVector& X, BlockedVector& Y, int level, bool InitialGuessIsZero) const
{
  if (level != NumLevels_ - 1)  // Perform one iteration of the V-cycle
  {
    // Apply presmoother
    SvecPre_[level]->Solve(X, Y, InitialGuessIsZero);

    // Compute residual
    BlockedVector DX = X.DeepCopy();
    Avec_[level]->Apply(Y, DX);
    DX.Update(1.0, X, -1.0);

    //  Create coarser representation of the residual
    int NV = X.GetVector(0)->NumVectors();
    Teuchos::RCP<BlockedVector> DXcoarse = Rvec_[level]->NewRangeBlockedVector(NV, false);
    Rvec_[level]->Apply(DX, *DXcoarse);

    // Damp error with coarser levels
    Teuchos::RCP<BlockedVector> DYcoarse = Pvec_[level]->NewDomainBlockedVector(NV, false);
    DoVcycle(*DXcoarse, *DYcoarse, level + 1, true);

    // Compute correction
    BlockedVector DY = Y.DeepCopy();
    Pvec_[level]->Apply(*DYcoarse, DY);
    Y.Update(1.0, DY, 1.0);

    // Apply post smoother
    SvecPos_[level]->Solve(X, Y, false);
  }
  else  // Apply presmoother
  {
    SvecPre_[level]->Solve(X, Y, InitialGuessIsZero);
  }


  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMGNXN::Vcycle::Solve(
    const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero) const
{
  // Check if everithing is set up
  if (!flag_set_up_A_) dserror("Operators missing");
  if (!flag_set_up_P_) dserror("Projectors missing");
  if (!flag_set_up_R_) dserror("Restrictors missing");
  if (!flag_set_up_Pre_) dserror("Pre-smoothers missing");
  if (!flag_set_up_Pos_) dserror("Post-smoothers missing");

  // Work!
  for (int i = 0; i < NumSweeps_; i++) DoVcycle(X, Y, FirstLevel_, InitialGuessIsZero and i == 0);
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMGNXN::VcycleSingle::VcycleSingle(int NumLevels, int NumSweeps, int FirstLevel)
    : NumLevels_(NumLevels),
      NumSweeps_(NumSweeps),
      FirstLevel_(FirstLevel),
      Avec_(NumLevels, Teuchos::null),
      Pvec_(NumLevels - 1, Teuchos::null),
      Rvec_(NumLevels - 1, Teuchos::null),
      SvecPre_(NumLevels, Teuchos::null),
      SvecPos_(NumLevels - 1, Teuchos::null),
      flag_set_up_A_(false),
      flag_set_up_P_(false),
      flag_set_up_R_(false),
      flag_set_up_Pre_(false),
      flag_set_up_Pos_(false)
{
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::VcycleSingle::SetOperators(
    std::vector<Teuchos::RCP<SparseMatrix>> Avec)
{
  if ((int)Avec.size() != NumLevels_) dserror("Error in Setting Avec_: Size dismatch.");
  for (int i = 0; i < NumLevels_; i++)
  {
    if (Avec[i] == Teuchos::null) dserror("Error in Setting Avec_: Null pointer.");
    Avec_[i] = Avec[i];
  }
  flag_set_up_A_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::VcycleSingle::SetProjectors(
    std::vector<Teuchos::RCP<SparseMatrix>> Pvec)
{
  if ((int)Pvec.size() != NumLevels_ - 1) dserror("Error in Setting Pvec_: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (Pvec[i] == Teuchos::null) dserror("Error in Setting Pvec_: Null pointer.");
    Pvec_[i] = Pvec[i];
  }
  flag_set_up_P_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::VcycleSingle::SetRestrictors(
    std::vector<Teuchos::RCP<SparseMatrix>> Rvec)
{
  if ((int)Rvec.size() != NumLevels_ - 1) dserror("Error in Setting Rvec_: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (Rvec[i] == Teuchos::null) dserror("Error in Setting Rvec_: Null pointer.");
    Rvec_[i] = Rvec[i];
  }
  flag_set_up_R_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::VcycleSingle::SetPreSmoothers(
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre)
{
  if ((int)SvecPre.size() != NumLevels_) dserror("Error in Setting SvecPre: Size dismatch.");
  for (int i = 0; i < NumLevels_; i++)
  {
    if (SvecPre[i] == Teuchos::null) dserror("Error in Setting SvecPre: Null pointer.");
    SvecPre_[i] = SvecPre[i];
  }
  flag_set_up_Pre_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::VcycleSingle::SetPosSmoothers(
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos)
{
  if ((int)SvecPos.size() != NumLevels_ - 1) dserror("Error in Setting SvecPos: Size dismatch.");
  for (int i = 0; i < NumLevels_ - 1; i++)
  {
    if (SvecPos[i] == Teuchos::null) dserror("Error in Setting SvecPos: Null pointer.");
    SvecPos_[i] = SvecPos[i];
  }
  flag_set_up_Pos_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMGNXN::VcycleSingle::DoVcycle(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level, bool InitialGuessIsZero) const
{
  if (level != NumLevels_ - 1)  // Perform one iteration of the V-cycle
  {
    // Apply presmoother
    SvecPre_[level]->Apply(X, Y, InitialGuessIsZero);

    // Compute residual
    int NV = X.NumVectors();
    Epetra_MultiVector DX(X.Map(), NV, false);
    Avec_[level]->Apply(Y, DX);
    DX.Update(1.0, X, -1.0);

    //  Create coarser representation of the residual
    const Epetra_Map& Map = Rvec_[level]->RangeMap();
    Epetra_MultiVector DXcoarse(Map, NV, false);
    Rvec_[level]->Apply(DX, DXcoarse);

    // Damp error with coarser levels
    const Epetra_Map& Map2 = Pvec_[level]->DomainMap();
    Epetra_MultiVector DYcoarse(Map2, NV, false);
    DoVcycle(DXcoarse, DYcoarse, level + 1, true);

    // Compute correction
    Epetra_MultiVector DY(Y.Map(), NV, false);
    Pvec_[level]->Apply(DYcoarse, DY);
    Y.Update(1.0, DY, 1.0);

    // Apply post smoother
    SvecPos_[level]->Apply(X, Y, false);
  }
  else  // Apply presmoother
  {
    SvecPre_[level]->Apply(X, Y, InitialGuessIsZero);
  }

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMGNXN::VcycleSingle::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  // Check if everithing is set up
  if (!flag_set_up_A_) dserror("Operators missing");
  if (!flag_set_up_P_) dserror("Projectors missing");
  if (!flag_set_up_R_) dserror("Restrictors missing");
  if (!flag_set_up_Pre_) dserror("Pre-smoothers missing");
  if (!flag_set_up_Pos_) dserror("Post-smoothers missing");

  // Work!
  for (int i = 0; i < NumSweeps_; i++) DoVcycle(X, Y, FirstLevel_, InitialGuessIsZero and i == 0);
  return;
}

#endif  // HAVE_MueLu
