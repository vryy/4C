/*----------------------------------------------------------------------*/
/*!
\file optimizer_graddesc.cpp

\brief Gradient descent with linesearch

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "optimizer_graddesc.H"
#include "invana_base.H"

#include "../linalg/linalg_utils.H"
#include "invana_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_comm/comm_utils.H"

/*----------------------------------------------------------------------*/
/* constructor                                                keh 01/13 */
/*----------------------------------------------------------------------*/
INVANA::OptimizerGradDesc::OptimizerGradDesc(const Teuchos::ParameterList& invp):
OptimizerBase(invp),
p_(Teuchos::null),
step_(Teuchos::null)
{;}

/*----------------------------------------------------------------------*/
/* setup algorithm specific stuff */
void INVANA::OptimizerGradDesc::Setup()
{
  p_= Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), 1, true));
  step_= Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), 1, true));
}


/*----------------------------------------------------------------------*/
/* optimization loop                                          keh 01/13 */
/*----------------------------------------------------------------------*/
void INVANA::OptimizerGradDesc::Integrate()
{
  if (!IsInit()) dserror("OptimizerBase is not inititialzed. Call Init() first");

  int success=0;

  // solve initially to get quantities:
  Evaluate(GetObjFunctVal(),GetGradient());

  //check gradient by fd:
#if 0
    std::cout << "gradient: " << *objgrad_ << std::endl;
    EvaluateGradientFD();
    std::cout << "gradient approx: " << *objgrad_ << std::endl;
#endif

  // old is current now
  UpdateGradient();

  //get search direction from gradient:
  p_->Update(-1.0, GetGradientView(), 0.0);

  UpdateObjFunctValue();

  MVNorm(GetGradientView(),SolLayoutMap(),2,&convcritc_);

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
    MVNorm(GetGradientView(),SolLayoutMap(),2,&convcritc_);

    //compute new direction only for runs
    p_->Update(-1.0, GetGradientView(), 0.0);

    // bring quantities to the next run
    UpdateGradient();
    UpdateObjFunctValue();
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
int INVANA::OptimizerGradDesc::EvaluateArmijoRule(double* tauopt, int* numsteps)
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

  MVNorm(GetGradientOldView(),SolLayoutMap(),2,&gnorm);

  double tau_n=std::min(1.0, 100/(1+gnorm));
  //std::cout << "trial step size: " << tau_n << std::endl;

  while (i<imax && tau_n<tau_max)
  {
    // step based on current stepsize
    step_->Update(tau_n, *p_, 0.0);

    //make a step
    UpdateSolution(*step_);

    Evaluate(GetObjFunctVal(),GetGradient());

    // check sufficient decrease:
    double dfp_o=0.0;
    MVDotProduct(GetGradientOldView(),*p_,SolLayoutMap(),&dfp_o);

    if ( (GetObjFunctValView()-GetObjFunctValOldView()) < c1*tau_n*dfp_o )
    {
      *tauopt=tau_n;
      *numsteps=i+1;
      return 0;
    }

    // do stepsize prediction based on polynomial models
    if (i==0)
      success=polymod(GetObjFunctValOldView(), dfp_o,tau_n,GetObjFunctValView(),blow,bhigh,tauopt);
    else
      success=polymod(GetObjFunctValOldView(),dfp_o,tau_n,GetObjFunctValView(),blow,bhigh,tau_l,e_l,tauopt);

    // repeat if cubic model fails
    if (success==1)
      success=polymod(GetObjFunctValOldView(), dfp_o,tau_n,GetObjFunctValView(),blow,bhigh,tauopt);

    e_l=GetObjFunctValView();
    tau_l=tau_n;
    tau_n=*tauopt;

    PrintLSStep(tau_l,i);

    UndoUpdateSolution();
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
int INVANA::OptimizerGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
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
int INVANA::OptimizerGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt)
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
void INVANA::OptimizerGradDesc::PrintOptStep(double tauopt, int numsteps)
{
  if (OptProb()->Comm().MyPID())
    return;

  printf("OPTIMIZATION STEP %3d | ", runc_);
  printf("Objective function: %10.8e | ", GetObjFunctValOldView());
  printf("Gradient : %10.8e | ", convcritc_);
  printf("stepsize : %10.8e | LSsteps %2d\n", tauopt, numsteps);
  fflush(stdout);

}

/*----------------------------------------------------------------------*/
/* print line search step */
void INVANA::OptimizerGradDesc::PrintLSStep(double tauopt, int numstep)
{

  if (OptProb()->Comm().MyPID()==0)
  {
    printf("   LINE SEARCH STEP %3d | ", numstep);
    printf("Objective function: %10.8e | ", GetObjFunctValView());
    printf("stepsize : %10.8e\n", tauopt);
    fflush(stdout);
  }

}

/*----------------------------------------------------------------------*/
/* print final results                                       kehl 01/13 */
/*----------------------------------------------------------------------*/
void INVANA::OptimizerGradDesc::Summarize()
{
  //std::cout << "the final vector of parameters: " << std::endl;
  //std::cout << *(Matman()->GetParams()) << std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/* Read restart                                               keh 03/14 */
/*----------------------------------------------------------------------*/
void INVANA::OptimizerGradDesc::ReadRestart(int run)
{
  IO::DiscretizationReader reader(OptProb()->Discret(),RestartFromFile(),run);

  if (not OptProb()->Comm().MyPID())
    std::cout << "Reading invana restart from step " << run << " from file: " << RestartFromFile()->FileName() << std::endl;

  runc_ = run;

  reader.ReadMultiVector(GetSolution(),"solution");

  return;
}


/*----------------------------------------------------------------------*/
/* Write restart                                              keh 03/14 */
/*----------------------------------------------------------------------*/
void INVANA::OptimizerGradDesc::WriteRestart()
{
  if (not OptProb()->Comm().MyPID())
    std::cout << "Writing invana restart for step " << runc_ <<  std::endl;

  Writer()->WriteNewStep(runc_);

  Writer()->WriteNamedVector("solution", GetSolution());

  return;
}
