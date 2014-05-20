/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_graddesc.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "../linalg/linalg_utils.H"
#include "invana_utils.H"
#include "stat_inv_ana_graddesc.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_adapter/ad_str_structure.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* constructor                                                keh 01/13 */
/*----------------------------------------------------------------------*/
STR::INVANA::StatInvAnaGradDesc::StatInvAnaGradDesc(Teuchos::RCP<DRT::Discretization> dis):
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

  p_= Teuchos::rcp(new Epetra_MultiVector(*(Matman()->ParamLayoutMap()), Matman()->NumVectors(),true));
  step_= Teuchos::rcp(new Epetra_MultiVector(*(Matman()->ParamLayoutMap()), Matman()->NumVectors(), true));

}


/*----------------------------------------------------------------------*/
/* optimization loop                                          keh 01/13 */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnaGradDesc::Optimize()
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

  // old is current now
  UpdateGradient();

  //get search direction from gradient:
  p_->Update(-1.0, *GetGradient(), 0.0);

  objval_o_ = objval_;

  MVNorm(GetGradient(),2,&convcritc_,Matman()->ParamLayoutMapUnique());

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
    MVNorm(GetGradient(),2,&convcritc_,Matman()->ParamLayoutMapUnique());

    //compute new direction only for runs
    p_->Update(-1.0, *GetGradient(), 0.0);

    // bring quantities to the next run
    UpdateGradient();
    objval_o_=objval_;
    runc_++;

    if (restartevry_ and (runc_%restartevry_ == 0))
        WriteRestart();

    //do some on screen printing
    PrintOptStep(tauopt, numsteps);
  }

  Summarize();

  return;
}


/*----------------------------------------------------------------------*/
/* do a line search based on armijo rule                      keh 03/14 */
/*----------------------------------------------------------------------*/
int STR::INVANA::StatInvAnaGradDesc::EvaluateArmijoRule(double* tauopt, int* numsteps)
{
  int i=0;
  int imax=20;
  double c1=1.0e-4;
  double tau_max=1.0e10;
  double gnorm=0.0;

  // "last"/"intermediate" values for cubic model
  // these are actually safeguardly set after the first call to the quadratic model
  double tau_l=0.0;
  double e_l=0.0;

  int success=0;

  //safeguard multiplicators
  double blow=0.1;
  double bhigh=0.5;

  MVNorm(GetGradientOld(),2,&gnorm,Matman()->ParamLayoutMapUnique());

  double tau_n=std::min(1.0, 100/(1+gnorm));
  //std::cout << "trial step size: " << tau_n << std::endl;

  while (i<imax && tau_n<tau_max)
  {
    // step based on current stepsize
    step_->Update(tau_n, *p_, 0.0);

    //make a step
    Matman()->UpdateParams(step_);
    SolveForwardProblem();
    SolveAdjointProblem();
    EvaluateGradient();
    EvaluateError();

    // check sufficient decrease:
    double dfp_o=0.0;
    MVDotProduct(GetGradientOld(),p_,&dfp_o,Matman()->ParamLayoutMapUnique());

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

    // repeat if cubic model fails
    if (success==1)
      success=polymod(objval_o_, dfp_o,tau_n,objval_,blow,bhigh,tauopt);

    e_l=objval_;
    tau_l=tau_n;
    tau_n=*tauopt;

    PrintLSStep(tau_l,i);

    Matman()->ResetParams();
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
/* quadratic model                                            keh 10/13 */
/*----------------------------------------------------------------------*/
int STR::INVANA::StatInvAnaGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  *tauopt=-dfp/(2*tau_n*(e_n-e_o-dfp));
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}


/*----------------------------------------------------------------------*/
/* cubic model                                               keh 10/13 */
/*----------------------------------------------------------------------*/
int STR::INVANA::StatInvAnaGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt)
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
/* print step information                                     keh 01/13 */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnaGradDesc::PrintOptStep(double tauopt, int numsteps)
{
  if (discret_->Comm().MyPID())
    return;

  printf("OPTIMIZATION STEP %3d | ", runc_);
  printf("Objective function: %10.8e | ", objval_o_);
  printf("Gradient : %10.8e | ", convcritc_);
  printf("stepsize : %10.8e | LSsteps %2d\n", tauopt, numsteps);
  fflush(stdout);

}

/*----------------------------------------------------------------------*/
/* print line search step */
void STR::INVANA::StatInvAnaGradDesc::PrintLSStep(double tauopt, int numstep)
{

  if (discret_->Comm().MyPID()==0)
  {
    printf("   LINE SEARCH STEP %3d | ", numstep);
    printf("Objective function: %10.8e | ", objval_);
    printf("stepsize : %10.8e\n", tauopt);
    fflush(stdout);
  }

}

/*----------------------------------------------------------------------*/
/* print final results                                       kehl 01/13 */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnaGradDesc::Summarize()
{
  if (not discret_->Comm().MyPID())
    std::cout << "the final vector of parameters: " << std::endl;

  std::cout << *(Matman()->GetParams()) << std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/* Read restart                                               keh 03/14 */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnaGradDesc::ReadRestart(int run)
{
  IO::DiscretizationReader reader(discret_,RestartFromFile(),run);

  if (not discret_->Comm().MyPID())
    std::cout << "Reading invana restart from step " << run << " from file: " << RestartFromFile()->FileName() << std::endl;

  //IO::DiscretizationReader reader(discret_, RestartFromFile(),run);
  if (run != reader.ReadInt("run"))
    dserror("Optimization step on file not equal to given step");

  runc_ = run;

  Teuchos::RCP<Epetra_MultiVector> params = Teuchos::rcp(new Epetra_MultiVector(*(Matman()->GetParams())));
  reader.ReadMultiVector(params,"optimization_parameters");
  Matman()->ReplaceParams(params);

  return;
}


/*----------------------------------------------------------------------*/
/* Write restart                                              keh 03/14 */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnaGradDesc::WriteRestart()
{
  if (not discret_->Comm().MyPID())
    std::cout << "Writing invana restart for step " << runc_ <<  std::endl;

  Writer()->NewStep(runc_, double(runc_));
  Writer()->WriteInt("run", runc_);

  // write vectors with unique gids only
  Teuchos::RCP<Epetra_MultiVector> uniqueparams = Teuchos::rcp(new Epetra_MultiVector(*Matman()->ParamLayoutMapUnique(), Matman()->NumVectors(),false));
  LINALG::Export(*(Matman()->GetParams()), *uniqueparams);

  Writer()->WriteVector("optimization_parameters", uniqueparams);

  return;
}
