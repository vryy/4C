/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_lbfgs.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_ana_lbfgs.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/drt_validparameters.H"

// needed to deal with materials
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* constructor */
STR::INVANA::StatInvAnaLBFGS::StatInvAnaLBFGS(Teuchos::RCP<DRT::Discretization> dis):
  StatInvAnalysis(dis),
stepsize_(0.0),
maxiter_(0),
runc_(0),
convcritc_(0)
{
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // max number of iterations
  maxiter_ = invap.get<int>("MAXITER");

  //set stepsize for gradient scheme
  stepsize_ = invap.get<double>("STEPSIZE");

  //get tolerance
  convtol_ = invap.get<double>("CONVTOL");

  //initialize storage vectors
  ssize_=invap.get<int>("SIZESTORAGE");
  ssize_=ssize_*matman_->NumParams();
  actsize_=0;

  sstore_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(-ssize_+1, 0, discret_->ElementColMap(), true));
  ystore_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(-ssize_+1, 0, discret_->ElementColMap(), true));

  p_= Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()), matman_->NumParams(),true));
  step_= Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()), matman_->NumParams(), true));

}

/*----------------------------------------------------------------------*/
/* do the optimization loop*/
void STR::INVANA::StatInvAnaLBFGS::Optimize()
{
  int success=0;

  // solve initially to get quantities:
  SolveForwardProblem();
  SolveAdjointProblem();
  EvaluateGradient();
  EvaluateError();

  //check gradient by fd:
#if 0
    std::cout << "gradient: " << *objgrad_ << std::endl;
    EvaluateGradientFD();
    std::cout << "gradient approx: " << *objgrad_ << std::endl;
#endif

  objgrad_o_->Update(1.0, *objgrad_, 0.0);

  //get search direction from gradient:
  p_->Update(-1.0, *objgrad_o_, 0.0);

  //objval_o_ = objfunct_->Evaluate(dis_,matman_);
  objval_o_ = objval_;


  MVNorm(objgrad_o_,2,&convcritc_);

  PrintOptStep(0,0);

  while (convcritc_ > convtol_ && runc_<maxiter_)
  {
    double tauopt=1.0;
    int numsteps=0;

    // do the line search
    success = EvaluateArmijoRule(&tauopt, &numsteps);

    if (success == 1)
    {
      std::cout << " Line Search Break Down" << std::endl;
      break;
    }

    //get the L2-norm:
    MVNorm(objgrad_,2,&convcritc_);

    //compute new direction only for runs>0
    if (runc_<1)
    {
      p_->Update(-1.0, *objgrad_, 0.0);
    }
    else
    {
      StoreVectors();
      ComputeDirection();
    }

    // bring quantities to the next run
    objgrad_o_->Update(1.0, *objgrad_, 0.0);
    objval_o_=objval_;
    runc_++;

    //do some on screen printing
    PrintOptStep(tauopt, numsteps);
  }

  Summarize();

  return;
}

/*----------------------------------------------------------------------*/
/* do a line search based on armijo rule */
int STR::INVANA::StatInvAnaLBFGS::EvaluateArmijoRule(double* tauopt, int* numsteps)
{
  int i=0;
  int imax=20;
  double c1=1.0e-4;
  double tau_max=1.0e10;
  double gnorm=0.0;

  // "last"/"intermediate" values for cubic model:
  double tau_l;
  double e_l;

  int success=0;

  //safeguard multiplicators
  double blow=0.1;
  double bhigh=0.5;

  MVNorm(objgrad_o_,2,&gnorm);

  double tau_n=std::min(1.0, 100/(1+gnorm));
  //std::cout << "trial step size: " << tau_n << std::endl;

  while (i<imax && tau_n<tau_max)
  {
    // step based on current stepsize
    step_->Update(tau_n, *p_, 0.0);

    //make a step
    matman_->UpdateParams(step_);
    SolveForwardProblem();
    SolveAdjointProblem();
    EvaluateGradient();
    //objval_ = objfunct_->Evaluate(dis_,matman_);
    EvaluateError();

    // check sufficient decrease:
    double dfp_o=0.0;
    MVDotProduct(objgrad_o_,p_,&dfp_o);

    if ( (objval_-objval_o_) < c1*tau_n*dfp_o )
    {
      *tauopt=tau_n;
      *numsteps=i+1;
      return 0;
    }

    // do stepsize prediction based on polynomial models
    if (i==0)
      success=polymod(objval_o_, dfp_o,tau_n,objval_,blow,bhigh,tauopt);
    else
      success=polymod(objval_o_,dfp_o,tau_n,objval_,blow,bhigh,tau_l,e_l,tauopt);

    if (success==1) return 1;

    e_l=objval_;
    tau_l=tau_n;
    tau_n=*tauopt;
    matman_->ResetParams();
    i++;

#if 0
    //brute force sampling:
    if (runc_==-1)
    {
      for (int k=0; k<100; k++)
      {
        step_->Update(tau_n*k/100, *p_, 0.0);
        matman_->UpdateParams(step_);
        SolveForwardProblem();
        double ee = objfunct_->Evaluate(dis_,matman_);
        std::cout << " run " << tau_n*k/100 << " " << ee << endl;
        matman_->ResetParams();
      }
    }
#endif
  }

  return 1;
}

/*----------------------------------------------------------------------*/
/* quadratic model */
int STR::INVANA::StatInvAnaLBFGS::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  *tauopt=-dfp/(2*tau_n*(e_n-e_o-dfp));
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
/* cubic model model */
int STR::INVANA::StatInvAnaLBFGS::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  double a1=tau_n*tau_n;
  double a2=tau_n*tau_n*tau_n;
  double a3=tau_l*tau_l;
  double a4=tau_l*tau_l*tau_l;

  double deta=a1*a4-a2*a3;

  if (deta<1.0e-14)
    return 1;

  double b1=e_n-(e_o+dfp*tau_n);
  double b2=e_l-(e_o+dfp*tau_l);

  double c1=1/deta*(a4*b1-a2*b2);
  double c2=1/deta*(-a3*b1+a1*b2);

  *tauopt=(-c1+sqrt(c1*c1-3*c2*dfp))/(3*c2);
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}


/*----------------------------------------------------------------------*/
/* store vectors*/
void STR::INVANA::StatInvAnaLBFGS::StoreVectors()
{
  if (runc_*matman_->NumParams()<=ssize_) // we have "<=" since we do not store the first run
    actsize_+=matman_->NumParams();

  Epetra_MultiVector s(*(discret_->ElementColMap()), (matman_->NumParams()),true);

  //push back s
  s.Update(1.0,*(matman_->GetParams()),-1.0,*(matman_->GetParamsOld()),0.0);
  for (int i=0; i<s.NumVectors(); i++)
    sstore_->UpdateSteps(*s(i));

  // push back y
  s.Update(1.0,*objgrad_,-1.0,*objgrad_o_,0.0);
  for (int i=0; i<s.NumVectors(); i++)
    ystore_->UpdateSteps(*s(i));


  return;
}

/*----------------------------------------------------------------------*/
/* compute new direction*/
void STR::INVANA::StatInvAnaLBFGS::ComputeDirection()
{
  p_->Update(1.0,*objgrad_,0.0);
  std::vector<double> alpha;

  // loop steps
  for (int i=0; i>-actsize_; i-=matman_->NumParams())
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;

    int ind=0;
    for (int j=matman_->NumParams(); j>0; j--)
    {
      //(*ystore_)(i-j+1)->Dot(*(*sstore_)(i-j+1),&a);
      //(*sstore_)(i-j+1)->Dot(*(*p_)(ind),&b);
      //MVDotProduct(Teuchos::rcp_dynamic_cast<Epetra_MultiVector>((*ystore_)(i-j+1)),Teuchos::rcp_dynamic_cast<Epetra_MultiVector>((*sstore_)(i-j+1)),&a);
      //MVDotProduct(Teuchos::rcp_dynamic_cast<Epetra_MultiVector>((*sstore_)(i-j+1)),Teuchos::rcp_dynamic_cast<Epetra_MultiVector>(Teuchos::rcp((*p_)(ind))),&b);
      MVDotProduct((*ystore_)(i-j+1),(*sstore_)(i-j+1),&a);
      MVDotProduct((*sstore_)(i-j+1),Teuchos::rcp((*p_)(ind), false),&b);
      ind++;
      aa += a;
      bb += b;
    }
    alpha.push_back(1/aa*bb);

    ind=0;
    for (int j=matman_->NumParams(); j>0; j--)
    {
      (*p_)(ind)->Update(-1.0*alpha.back(), *(*ystore_)(i-j+1),1.0 );
      ind++;
    }
  }

  // Some saling of the initial hessian might come in here but has not been proven to be effective
  // altough they say so

  for (int i=-actsize_+1; i<=0; i+=matman_->NumParams())
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;
    double beta=0.0;

    for (int j=0; j<matman_->NumParams(); j++)
    {
      //(*ystore_)(i+j)->Dot(*(*sstore_)(i+j),&a);
      //(*ystore_)(i+j)->Dot(*(*p_)(j),&b);
      MVDotProduct((*ystore_)(i+j),(*sstore_)(i+j),&a);
      MVDotProduct((*ystore_)(i+j),Teuchos::rcp((*p_)(j), false),&b);
      aa += a;
      bb += b;
    }

    beta=1/aa*bb;
    double alphac=alpha.back();
    alpha.pop_back();

    for (int j=0; j<matman_->NumParams(); j++)
      (*p_)(j)->Update(alphac-beta, *(*sstore_)(i+j),1.0 );
  }

  // we do minimization not maximization
  p_->Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* print final results*/
void STR::INVANA::StatInvAnaLBFGS::PrintOptStep(double tauopt, int numsteps)
{
  if (discret_->Comm().MyPID()==0)
  {
    printf("OPTIMIZATION STEP %3d | ", runc_);
    printf("Objective function: %10.8e | ", objval_o_);
    printf("Gradient : %10.8e | ", convcritc_);
    printf("stepsize : %10.8e | LSsteps %2d\n", tauopt, numsteps);
    fflush(stdout);
  }

}

/*----------------------------------------------------------------------*/
/* print final results*/
void STR::INVANA::StatInvAnaLBFGS::Summarize()
{
  std::cout << "the final vector of parameters: " << std::endl;
  std::cout << *(matman_->GetParams()) << std::endl;
  return;
}

