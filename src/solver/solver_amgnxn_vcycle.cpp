/*!----------------------------------------------------------------------
\file solver_amgnxn_vcycle.cpp

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
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
LINALG::SOLVER::Richardson_Vcycle_Operator::Richardson_Vcycle_Operator
(int NumLevels,int NumSweeps,double omega):
NumLevels_(NumLevels),
NumSweeps_(NumSweeps),
omega_(omega),
Avec_(NumLevels,Teuchos::null),
Pvec_(NumLevels-1,Teuchos::null),
Rvec_(NumLevels-1,Teuchos::null),
SvecPre_(NumLevels,Teuchos::null),
SvecPos_(NumLevels-1,Teuchos::null),
flag_set_up_A_(false),
flag_set_up_P_(false),
flag_set_up_R_(false),
flag_set_up_Pre_(false),
flag_set_up_Pos_(false)
{}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Richardson_Vcycle_Operator::SetOperators
  (std::vector< Teuchos::RCP<Epetra_Operator> > Avec)
{
  if((int)Avec.size()!=NumLevels_)
    dserror("Error in Setting Avec_: Size dismatch.");
  for(int i=0;i<NumLevels_;i++)
  {
    if(Avec[i]==Teuchos::null)
      dserror("Error in Setting Avec_: Null pointer.");
    Avec_[i]=Avec[i];
  }
  flag_set_up_A_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Richardson_Vcycle_Operator::SetProjectors
  (std::vector< Teuchos::RCP<Epetra_Operator> > Pvec)
{
  if((int)Pvec.size()!=NumLevels_-1)
    dserror("Error in Setting Pvec_: Size dismatch.");
  for(int i=0;i<NumLevels_-1;i++)
  {
    if(Pvec[i]==Teuchos::null)
      dserror("Error in Setting Pvec_: Null pointer.");
    Pvec_[i]=Pvec[i];
  }
  flag_set_up_P_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Richardson_Vcycle_Operator::SetRestrictors
  (std::vector< Teuchos::RCP<Epetra_Operator> > Rvec)
{
  if((int)Rvec.size()!=NumLevels_-1)
    dserror("Error in Setting Rvec_: Size dismatch.");
  for(int i=0;i<NumLevels_-1;i++)
  {
    if(Rvec[i]==Teuchos::null)
      dserror("Error in Setting Rvec_: Null pointer.");
    Rvec_[i]=Rvec[i];
  }
  flag_set_up_R_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Richardson_Vcycle_Operator::SetPreSmoothers
  (std::vector< Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase> > SvecPre)
{
  if((int)SvecPre.size()!=NumLevels_)
    dserror("Error in Setting SvecPre: Size dismatch.");
  for(int i=0;i<NumLevels_;i++)
  {
    if(SvecPre[i]==Teuchos::null)
      dserror("Error in Setting SvecPre: Null pointer.");
    SvecPre_[i]=SvecPre[i];
  }
  flag_set_up_Pre_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Richardson_Vcycle_Operator::SetPosSmoothers
  (std::vector< Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase> > SvecPos)
{
  if((int)SvecPos.size()!=NumLevels_-1)
    dserror("Error in Setting SvecPos: Size dismatch.");
  for(int i=0;i<NumLevels_-1;i++)
  {
    if(SvecPos[i]==Teuchos::null)
      dserror("Error in Setting SvecPos: Null pointer.");
    SvecPos_[i]=SvecPos[i];
  }
  flag_set_up_Pos_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::Richardson_Vcycle_Operator::Vcycle(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level, bool InitialGuessIsZero) const
{


  if(level!=NumLevels_-1) // Perform one iteration of the V-cycle
  {
    // Apply presmoother
    SvecPre_[level]->Apply(X,Y,InitialGuessIsZero);

    // Compute residual TODO optimize if InitialGuessIsZero == true
    Epetra_MultiVector DX(X.Map(),X.NumVectors());
    Avec_[level]->Apply(Y,DX);
    DX.Update(1.0,X,-1.0);

    //  Create coarser representation of the residual
    Epetra_MultiVector DXcoarse(Rvec_[level]->OperatorRangeMap(),X.NumVectors());
    Rvec_[level]->Apply(DX,DXcoarse);

    // Damp error with coarser levels
    Epetra_MultiVector DYcoarse(Pvec_[level]->OperatorDomainMap(),X.NumVectors(),true);
    DYcoarse.PutScalar(0.0);
    Vcycle(DXcoarse,DYcoarse,level+1,true);

    // Compute correction
    Epetra_MultiVector DY(Y.Map(),X.NumVectors());
    Pvec_[level]->Apply(DYcoarse,DY);
    Y.Update(1.0,DY,1.0);

    // Apply post smoother
    SvecPos_[level]->Apply(X,Y,false);

  }
  else // Apply presmoother
  {
    Epetra_MultiVector X00(X.Map(),X.NumVectors());
    SvecPre_[level]->Apply(X,Y,InitialGuessIsZero);
  }

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::Richardson_Vcycle_Operator::Richardson_Vcycle(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, int start_level) const
{
  // Create auxiliary vectors
  Epetra_MultiVector Ytmp(Y.Map(),X.NumVectors(),true); // true is necessary! (initial guess has to be zero)
  Epetra_MultiVector DX(X.Map(),X.NumVectors(),false);
  Epetra_MultiVector DY(Y.Map(),X.NumVectors(),false);

  for(int i=0;i<NumSweeps_;i++)
  {

    double scal_aux = (i==0)? 0.0 : 1.0;

    // Compute residual
    if(i!=0)
      Avec_[0]->Apply(Ytmp,DX);
    DX.Update(1.0,X,-1.0*scal_aux);

    // Apply V-cycle as preconditioner
    DY.PutScalar(0.0);
    Vcycle(DX,DY,start_level,true);

    // Apply correction
    Ytmp.Update(omega_,DY,scal_aux);

  }
  Y=Ytmp;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::Richardson_Vcycle_Operator::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, int start_level) const
{
  // Check if everithing is set up
  if(!flag_set_up_A_)
    dserror("Operators missing");
  if(!flag_set_up_P_)
    dserror("Projectors missing");
  if(!flag_set_up_R_)
    dserror("Restrictors missing");
  if(!flag_set_up_Pre_)
    dserror("Pre-smoothers missing");
  if(!flag_set_up_Pos_)
    dserror("Post-smoothers missing");

  // Work!
  Richardson_Vcycle(X,Y,start_level);
  return;
}



#endif // HAVE_MueLu
