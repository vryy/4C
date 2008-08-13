/*!----------------------------------------------------------------------
\file uzawasolver.cpp

\brief Class constaining uzawa algorithm to solve linear system.

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "uzawasolver.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
UTILS::UzawaSolver::UzawaSolver(
    RCP<DRT::Discretization> discr,
    LINALG::Solver& solver,
    RCP<Epetra_Vector>    dirichtoggle,
    RCP<Epetra_Vector>    invtoggle,
    ParameterList params):
actdisc_(discr),
dirichtoggle_(&(*dirichtoggle),false),
invtoggle_(&(*invtoggle),false)
{
  solver_         = rcp(&solver,false);
  toldisp_        = params.get<double>("tolerance displacements", 1.0e-07);
  tolres_         = params.get<double>("tolerance residual", 1.0e-07);
  tolconstr_      = params.get<double>("tolerance constraint"     ,1.0e-07);
  // this exception handler is not the nicest thing, 
  // but the original (and verified) input 
  // parameter list is copied to another parameter list in case of 
  // the StruGenAlpha time integrator object (and apparently other time 
  // integrators). This is actually done in dyn_nlnstructural_drt(). 
  // The approach in StruTimIntImpl (and thus StruTimIntGenAlpha)
  // is not to change the parameter list after input nor to copy
  // them to new lists. This copy mechanism does not improve the
  // data. 
  // Thus we need this exception handler to getting along.
  try
  {
    // for StruGenAlpha
    isadapttol_   = params.get<bool>("ADAPTCONV",true);
  }
  catch (const Teuchos::Exceptions::InvalidParameterType)
  {
    // for StruTimIntImpl
    isadapttol_   = true;
    isadapttol_   = (Teuchos::getIntegralValue<int>(params,"ADAPTCONV") == 1);
  }
  adaptolbetter_  = params.get<double>("ADAPTCONV_BETTER",0.01);
  maxIter_        = params.get<int>   ("uzawa maxiter"         ,50);
  uzawaParam_     = params.get<double>("uzawa parameter",1);
  counter_        = 0;
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 11/07|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
*-----------------------------------------------------------------------*/
void UTILS::UzawaSolver::Solve(
        RCP<LINALG::SparseMatrix> stiff,
        RCP<LINALG::SparseMatrix> constr,
        RCP<Epetra_Vector> dispinc,
        RCP<Epetra_Vector> lagrinc,
        RCP<Epetra_Vector> rhsstand,
        RCP<Epetra_Vector> rhsconstr
        )
{
  
  // For every iteration step an uzawa algorithm is used to solve the linear system.
  //Preparation of uzawa method to solve the linear system.
  double norm_uzawa;
  double norm_uzawa_old;
  double quotient;
  double norm_constr_uzawa;
  int numiter_uzawa=0;
  //counter used for adaptivity
  int count_paramadapt=1;

  RCP<Epetra_Vector> constrTLagrInc=rcp(new Epetra_Vector(rhsstand->Map()));
  RCP<Epetra_Vector> constrTDispInc= rcp(new Epetra_Vector(rhsconstr->Map()));

  // Compute residual of the uzawa algorithm

  constr->Multiply(false,*lagrinc,*constrTLagrInc);//necessary?
  
  RCP<Epetra_Vector> fresmcopy=rcp(new Epetra_Vector(*rhsstand));
  fresmcopy->Update(1.0,*constrTLagrInc,1.0);
  Epetra_Vector uzawa_res(*fresmcopy);
  (*stiff).Multiply(false,*dispinc,uzawa_res);
  uzawa_res.Update(1.0,*fresmcopy,-1.0);
  
  // blank residual DOFs which are on Dirichlet BC 
  Epetra_Vector rescopy(uzawa_res);
  uzawa_res.Multiply(1.0,*invtoggle_,rescopy,0.0);
  
  uzawa_res.Norm2(&norm_uzawa);
  Epetra_Vector constr_res(lagrinc->Map());
  
  constr->Multiply(true,*dispinc,*constrTDispInc) ;
  constr_res.Update(1.0,*constrTDispInc,1.0,*(rhsconstr),0.0);
  constr_res.Norm2(&norm_constr_uzawa);
  quotient =1;
  RCP<Epetra_Vector> zeros = rcp(new Epetra_Vector(rhsstand->Map(),true));
  //Solve one iteration step with augmented lagrange
  //Since we calculate displacement norm as well, at least one step has to be taken
  while (((norm_uzawa > tolres_/10 or norm_constr_uzawa>tolconstr_/10)
          and numiter_uzawa < maxIter_) or numiter_uzawa<1)
  {
      LINALG::ApplyDirichlettoSystem(stiff,dispinc,fresmcopy,zeros,dirichtoggle_);
      // solve for disi
      // Solve K . IncD = -R  ===>  IncD_{n+1}
      if (isadapttol_ && counter_ && numiter_uzawa)
      {
        double worst = norm_uzawa;
        double wanted = tolres_/10.0;
        solver_->AdaptTolerance(wanted,worst,adaptolbetter_);
      }
      solver_->Solve(stiff->EpetraMatrix(),dispinc,fresmcopy,true,numiter_uzawa==0 && counter_==0);
      solver_->ResetTolerance();

      //compute Lagrange multiplier increment
      constrTDispInc->PutScalar(0.0);
      constr->Multiply(true,*dispinc,*constrTDispInc) ;
      lagrinc->Update(uzawaParam_,*constrTDispInc,uzawaParam_,*rhsconstr,1.0);
      
      //Compute residual of the uzawa algorithm
      constr->Multiply(false,*lagrinc,*constrTLagrInc);
      constrTLagrInc->Scale(-1.0);
      
      fresmcopy->Update(1.0,*constrTLagrInc,1.0,*rhsstand,0.0);
      Epetra_Vector uzawa_res(*fresmcopy);
      (*stiff).Multiply(false,*dispinc,uzawa_res);
      uzawa_res.Update(1.0,*fresmcopy,-1.0);
      // blank residual DOFs which are on Dirichlet BC
      Epetra_Vector rescopy(uzawa_res);
      uzawa_res.Multiply(1.0,*invtoggle_,rescopy,0.0);
      norm_uzawa_old=norm_uzawa;
      uzawa_res.Norm2(&norm_uzawa);
      Epetra_Vector constr_res(lagrinc->Map());

      constr_res.Update(1.0,*constrTDispInc,1.0,*rhsconstr,0.0);
      constr_res.Norm2(& norm_constr_uzawa);

      //-------------Adapt Uzawa parameter--------------
      // For a constant parameter the quotient of two successive residual norms
      // stays constant during the computation. So this quotient seems to be a good
      // measure for the parameter choice

      // Adaptivity only takes place every second step. Otherwise the quotient is not significant.
      if (count_paramadapt>=2)
      {
        double quotient_new=norm_uzawa/norm_uzawa_old;
        // In case of divergence the parameter must be too high
        if (quotient_new>1)
        {
          uzawaParam_=uzawaParam_/2.;
          count_paramadapt=0;
          quotient=1;
        }
        else
        {
          // In case the newly computed quotient is better than the one obtained from the
          // previous parameter, the parameter is increased by a factor (1+quotient_new)
          if (quotient>quotient_new)
          {
            uzawaParam_=uzawaParam_*(1+quotient_new);
            quotient=quotient_new;
            count_paramadapt=0;
          }
          // In case the newly computed quotient is worse than the one obtained from the
          // previous parameter, the parameter is decreased by a factor 1/(1+quotient_new)
          else
          {
            uzawaParam_=uzawaParam_/(1+quotient_new);
            quotient=quotient_new;
            count_paramadapt=0;
          }
        }
      }
      if (uzawaParam_<1.0E-8)
      {
        uzawaParam_*=1.0E3;
        break;
      }
      count_paramadapt++;
      numiter_uzawa++;
  } //Uzawa loop
  
  int myrank=(actdisc_->Comm().MyPID());
  if (!myrank)
  {
      cout<<"Uzawa steps "<<numiter_uzawa<<", Uzawa parameter: "<< uzawaParam_<< endl;
  }
  counter_++;
  return;
}


#endif
