/*!----------------------------------------------------------------------
\file constraintsolver.cpp

\brief Class containing uzawa algorithm to solve linear system.

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraintsolver.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
UTILS::ConstraintSolver::ConstraintSolver
(
  RCP<DRT::Discretization> discr,
  LINALG::Solver& solver,
  RCP<Epetra_Vector>    dirichtoggle,
  RCP<Epetra_Vector>    invtoggle,
  ParameterList params
):
actdisc_(discr),
maxIter_(params.get<int>   ("UZAWAMAXITER", 50)),
dirichtoggle_(&(*dirichtoggle),false),
invtoggle_(&(*invtoggle),false)
{
  solver_ = rcp(&solver,false);
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
  
  int algochoice;
  try
  {
    // for StruGenAlpha
    isadapttol_ = params.get<bool>("ADAPTCONV",true);
    algochoice = params.get<int>("UZAWAALGO");
  }
  catch (const Teuchos::Exceptions::InvalidParameterType)
  {
    // for StruTimIntImpl
    isadapttol_ = true;
    isadapttol_ = (Teuchos::getIntegralValue<int>(params,"ADAPTCONV") == 1);
    algochoice = getIntegralValue<int>(params,"UZAWAALGO");
  }
  adaptolbetter_ = params.get<double>("ADAPTCONV_BETTER", 0.01);
  uzawaparam_ = params.get<double>("UZAWAPARAM", 1);
  minparam_ = uzawaparam_*1E-3;
  uzawatol_ = params.get<double>("UZAWATOL", 1E-8);
  
  switch (algochoice)
  {
  case (INPUTPARAMS::consolve_iterative):
    algo_ = UTILS::ConstraintSolver::iterative;
  break;
  case (INPUTPARAMS::consolve_direct):
    algo_ = UTILS::ConstraintSolver::direct;
//    #ifdef PARALLEL
//      dserror("Direct constraint solver is not working in parallel!");
//    #endif
  break;
  default:
    dserror("Unknown type of constraint solver algorithm. Can be 'iterative' or 'direct'!");
  }
  counter_ = 0;
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system                                        |
*-----------------------------------------------------------------------*/
void UTILS::ConstraintSolver::Solve
(
  RCP<LINALG::SparseMatrix> stiff,
  RCP<LINALG::SparseMatrix> constr,
  RCP<Epetra_Vector> dispinc,
  RCP<Epetra_Vector> lagrinc,
  const RCP<Epetra_Vector> rhsstand,
  const RCP<Epetra_Vector> rhsconstr
)
{
  switch (algo_)
  {
    case iterative:
      SolveIterative(stiff,constr,dispinc,lagrinc,rhsstand,rhsconstr);
      break;
    case direct:
      SolveDirect(stiff,constr,dispinc,lagrinc,rhsstand,rhsconstr);
      break;
    default :
      dserror("Unknown constraint solution technique!");
  }
  return;
}
  
/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system by iterative Uzawa algorithm           |
*-----------------------------------------------------------------------*/
void UTILS::ConstraintSolver::SolveIterative
(
  RCP<LINALG::SparseMatrix> stiff,
  RCP<LINALG::SparseMatrix> constr,
  RCP<Epetra_Vector> dispinc,
  RCP<Epetra_Vector> lagrinc,
  const RCP<Epetra_Vector> rhsstand,
  const RCP<Epetra_Vector> rhsconstr
)
{ 
  const int myrank=(actdisc_->Comm().MyPID());
  // For every iteration step an uzawa algorithm is used to solve the linear system.
  //Preparation of uzawa method to solve the linear system.
  double norm_uzawa;
  double norm_uzawa_old;
  double quotient;
  double norm_constr_uzawa;
  int numiter_uzawa = 0;
  //counter used for adaptivity
  const int adaptstep = 2;
  const int minstep = 2;
  int count_paramadapt = 1;

  const double computol = 1E-8;
  
  RCP<Epetra_Vector> constrTLagrInc = rcp(new Epetra_Vector(rhsstand->Map()));
  RCP<Epetra_Vector> constrTDispInc = rcp(new Epetra_Vector(rhsconstr->Map()));
  
  RCP<Epetra_Vector> zeros = rcp(new Epetra_Vector(rhsstand->Map(),true));
  
  // Compute residual of the uzawa algorithm
  RCP<Epetra_Vector> fresmcopy=rcp(new Epetra_Vector(*rhsstand));
  Epetra_Vector uzawa_res(*fresmcopy);
  (*stiff).Multiply(false,*dispinc,uzawa_res);
  uzawa_res.Update(1.0,*fresmcopy,-1.0);
  
  // blank residual DOFs which are on Dirichlet BC 
  Epetra_Vector rescopy(uzawa_res);
  uzawa_res.Multiply(1.0,*invtoggle_,rescopy,0.0);
  
  uzawa_res.Norm2(&norm_uzawa);
  Epetra_Vector constr_res(lagrinc->Map());
  
  constr_res.Update(1.0,*(rhsconstr),0.0);
  constr_res.Norm2(&norm_constr_uzawa);
  quotient =1;
  //Solve one iteration step with augmented lagrange
  //Since we calculate displacement norm as well, at least one step has to be taken
  while (((norm_uzawa > uzawatol_ or norm_constr_uzawa > uzawatol_) and numiter_uzawa < maxIter_) 
      or numiter_uzawa < minstep)
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
    lagrinc->Update(uzawaparam_,*constrTDispInc,uzawaparam_,*rhsconstr,1.0);
    
    //Compute residual of the uzawa algorithm
    constr->Multiply(false,*lagrinc,*constrTLagrInc);
    
    fresmcopy->Update(-1.0,*constrTLagrInc,1.0,*rhsstand,0.0);
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
    constr_res.Norm2(&norm_constr_uzawa);
    //-------------Adapt Uzawa parameter--------------
    // For a constant parameter the quotient of two successive residual norms
    // stays nearly constant during the computation. So this quotient seems to be a good
    // measure for the parameter choice
    // Adaptivity only takes place every second step. Otherwise the quotient is not significant.
    if (count_paramadapt>=adaptstep)
    {
      double quotient_new=norm_uzawa/norm_uzawa_old;
      // In case of divergence the parameter must be too high
      if (quotient_new>(1.+computol))
      {
        if (uzawaparam_>2.*minparam_)
          uzawaparam_ = uzawaparam_/2.;
        quotient=1;
      }
      else 
      {
        // In case the newly computed quotient is better than the one obtained from the
        // previous parameter, the parameter is increased by a factor (1+quotient_new)
        if (quotient>=quotient_new)
        {
          uzawaparam_=uzawaparam_*(1.+quotient_new);
          quotient=quotient_new;
        }
        // In case the newly computed quotient is worse than the one obtained from the
        // previous parameter, the parameter is decreased by a factor 1/(1+quotient_new)
        else 
        {
          if (uzawaparam_>2.*minparam_)
            uzawaparam_=uzawaparam_/(1.+quotient_new);
          quotient=quotient_new;
        }
      }
      
      if (uzawaparam_<=minparam_)
      {
        if (!myrank)
          cout<<"leaving uzawa loop since Uzawa parameter is too low"<<endl;
        uzawaparam_*=1E2;
        break;
      }
      count_paramadapt=0;
    }
    count_paramadapt++;
    numiter_uzawa++;
  } //Uzawa loop
  
  if (!myrank)
  {
     cout<<"Uzawa steps "<<numiter_uzawa<<", Uzawa parameter: "<< uzawaparam_;
     //cout<<endl;
     cout<<", residual norms for linear system: "<< norm_constr_uzawa<<" and "<<norm_uzawa<<endl;
  }
  counter_++;
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system by iterative Uzawa algorithm           |
*-----------------------------------------------------------------------*/
void UTILS::ConstraintSolver::SolveDirect
(
  RCP<LINALG::SparseMatrix> stiff,
  RCP<LINALG::SparseMatrix> constr,
  RCP<Epetra_Vector> dispinc,
  RCP<Epetra_Vector> lagrinc,
  const RCP<Epetra_Vector> rhsstand,
  const RCP<Epetra_Vector> rhsconstr
)
{
  // define maps of standard dofs and additional lagrange multipliers
  RCP<Epetra_Map> standrowmap = rcp(new Epetra_Map(stiff->RowMap()));
  RCP<Epetra_Map> conrowmap = rcp(new Epetra_Map(constr->DomainMap()));
  // merge maps to one large map
  RCP<Epetra_Map> mergedmap = LINALG::MergeMap(standrowmap,conrowmap,false);
  
  // define MapExtractor
  LINALG::MapExtractor mapext(*mergedmap,standrowmap,conrowmap);
  
  // initialize large Sparse Matrix and Epetra_Vectors
  RCP<LINALG::SparseMatrix> mergedmatrix = rcp(new LINALG::SparseMatrix(*mergedmap,mergedmap->NumMyElements()));
  RCP<Epetra_Vector> mergedrhs = rcp(new Epetra_Vector(*mergedmap));
  RCP<Epetra_Vector> mergedsol = rcp(new Epetra_Vector(*mergedmap));
  RCP<Epetra_Vector> mergeddtog = rcp(new Epetra_Vector(*mergedmap));
  RCP<Epetra_Vector> mergedzeros = rcp(new Epetra_Vector(*mergedmap));
  
  // fill merged matrix using Add
  mergedmatrix -> Add(*stiff,false,1.0,1.0);
  mergedmatrix -> Add(*constr,false,1.0,1.0);
  mergedmatrix -> Add(*constr,true,1.0,1.0);
  mergedmatrix -> Complete(*mergedmap,*mergedmap);

  // fill merged vectors using Export
  LINALG::Export(*rhsconstr,*mergedrhs);
  mergedrhs -> Scale(-1.0);
  LINALG::Export(*rhsstand,*mergedrhs);
  LINALG::Export(*dirichtoggle_,*mergeddtog);

  // apply dirichlet boundary conditions
  LINALG::ApplyDirichlettoSystem(mergedmatrix,mergedsol,mergedrhs,mergedzeros,mergeddtog);
  
  // solve
  solver_->Solve(mergedmatrix->EpetraMatrix(),mergedsol,mergedrhs,true,counter_==0);
  solver_->ResetTolerance();

  // store results in smaller vectors
  mapext.ExtractCondVector(mergedsol,dispinc);
  mapext.ExtractOtherVector(mergedsol,lagrinc);
  
  counter_++;
  return;
}

#endif
