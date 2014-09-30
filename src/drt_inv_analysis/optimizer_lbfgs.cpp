/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_lbfgs.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "optimizer_lbfgs.H"
#include "invana_base.H"

#include "../linalg/linalg_utils.H"
#include "invana_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_comm/comm_utils.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::INVANA::OptimizerLBFGS::OptimizerLBFGS(const Teuchos::ParameterList& invp):
  OptimizerBase(invp),
sstore_(Teuchos::null),
ystore_(Teuchos::null),
p_(Teuchos::null),
step_(Teuchos::null)
{
  //initialize storage vectors
  ssize_=invp.get<int>("SIZESTORAGE");

}

/*----------------------------------------------------------------------*/
/* setup algorithm specific stuff */
void STR::INVANA::OptimizerLBFGS::Setup()
{
  ssize_=ssize_*numvecs_;
  actsize_=0;

  sstore_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, &SolLayoutMap(), true));
  ystore_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, &SolLayoutMap(), true));

  p_= Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), numvecs_, true));
  step_= Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), numvecs_, true));
}

/*----------------------------------------------------------------------*/
/* do the optimization loop*/
void STR::INVANA::OptimizerLBFGS::Integrate()
{
  if (!IsInit()) dserror("OpimizerBase is not inititialzed. Call Init() first");

  int success=0;

  // solve initially to get quantities:
  Evaluate(GetObjFunctVal(),GetGradient());

  //check gradient by fd:
#if 0
    std::cout << "gradient: " << std::endl;
    PrintDataToScreen(*GetGradient());
    EvaluateGradientFD();

    std::cout << "gradient approx: " << std::endl;
    PrintDataToScreen(*GetGradient());
    exit(0);
#endif

  UpdateGradient();

  //get search direction from gradient:
  p_->Update(-1.0, GetGradientOldView(), 0.0);
  if (DRT::Problem::Instance()->Restart())
    ComputeDirection();

  UpdateObjFunctValue();
  error_incr_ = GetObjFunctValView();

  MVNorm(GetGradientOldView(),SolLayoutMapUnique(),2,&convcritc_);

  PrintOptStep(0,0);

  while (convcritc_ > convtol_ && runc_<maxiter_)
  {
    double tauopt=1.0;
    int numsteps=0;

    // do the line search
    success = EvaluateArmijoRule(&tauopt, &numsteps);

    if (success == 1)
    {
      std::cout << " Line Search Break Down with current stepsize " << tauopt << std::endl;
      break;
    }

    //get the L2-norm:
    MVNorm(GetGradientView(),SolLayoutMapUnique(),2,&convcritc_);

    //compute new direction only for runs>0
    if (runc_<1)
    {
      p_->Update(-1.0, GetGradientView(), 0.0);
    }
    else
    {
      StoreVectors();
      ComputeDirection();
    }

    // bring quantities to the next run
    UpdateGradient();
    error_incr_=GetObjFunctValOldView()-GetObjFunctValView();
    UpdateObjFunctValue();
    runc_++;

    //do some on screen printing
    PrintOptStep(tauopt, numsteps);

    if (restartevry_ and (runc_%restartevry_ == 0) and runc_)
      WriteRestart();
  }

  Evaluate(NULL,Teuchos::null); // set optimization parameters to the optproblem
  Summarize();

  return;
}

/*----------------------------------------------------------------------*/
/* do a line search based on armijo rule */
int STR::INVANA::OptimizerLBFGS::EvaluateArmijoRule(double* tauopt, int* numsteps)
{
  int i=0;
  int imax=20;
  double c1=1.0e-4;
  //double c2=0.4;
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

  MVNorm(GetGradientOldView(),SolLayoutMapUnique(),2,&gnorm);

  double tau_n=std::min(stepsize_, 100.0/(1+gnorm));
  //std::cout << "trial step size: " << tau_n << std::endl;

  // step based on current stepsize
  step_->Update(tau_n, *p_, 0.0);

  while (i<imax && tau_n<tau_max)
  {
    //make a step
    UpdateSolution(*step_);

    Evaluate(GetObjFunctVal(),GetGradient());

    // check sufficient decrease:
    double dfp_o=0.0;
    MVDotProduct(GetGradientOldView(),*p_,SolLayoutMapUnique(),&dfp_o);

    if ( (GetObjFunctValView()-GetObjFunctValOldView()) < c1*tau_n*dfp_o )
    {
      double dfp_n=0.0;
      MVDotProduct(GetGradientView(),*p_,SolLayoutMapUnique(),&dfp_n);
      //check curvature condition:
      if (1) //(abs(dfp_n)<c2*abs(dfp_o))
      {
        *tauopt=tau_n;
        *numsteps=i+1;
        return 0;
      }
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

    // next step based on proposed tauopt;
    step_->Update(tau_n, *p_, 0.0);

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
        Matman()->UpdateParams(step_);
        SolveForwardProblem();
        double ee = objfunct_->Evaluate(dis_,Matman());
        std::cout << " run " << tau_n*k/100 << " " << ee << endl;
        Matman()->ResetParams();
      }
    }
#endif
  }

  return 1;
}

/*----------------------------------------------------------------------*/
/* quadratic model */
int STR::INVANA::OptimizerLBFGS::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  *tauopt=-(dfp*tau_n*tau_n)/(2*(e_n-e_o-dfp*tau_n));
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
/* cubic model model */
int STR::INVANA::OptimizerLBFGS::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  double a1=tau_n*tau_n;
  double a2=tau_n*tau_n*tau_n;
  double a3=tau_l*tau_l;
  double a4=tau_l*tau_l*tau_l;

  double deta=a1*a4-a2*a3;

  if (deta<1.0e-14) return 1;

  double b1=e_n-(e_o+dfp*tau_n);
  double b2=e_l-(e_o+dfp*tau_l);

  double c1=1/deta*(a4*b1-a2*b2);
  double c2=1/deta*(-a3*b1+a1*b2);

  double arg=c1*c1-3*c2*dfp;

  if (arg < 0.0) return 1;

  *tauopt=(-c1+sqrt(arg))/(3*c2);
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
/* store vectors*/
void STR::INVANA::OptimizerLBFGS::StoreVectors()
{
  if (runc_*numvecs_<=ssize_) // we have "<=" since we do not store the first run
    actsize_+=numvecs_;

  Epetra_MultiVector s(SolLayoutMap(), (numvecs_),true);

  //push back s
  s.Update(1.0,GetSolutionView(),-1.0,GetSolutionOldView(),0.0);
  for (int i=0; i<s.NumVectors(); i++)
    sstore_->UpdateSteps(*s(i));

  // push back y
  s.Update(1.0,GetGradientView(),-1.0,GetGradientOldView(),0.0);
  for (int i=0; i<s.NumVectors(); i++)
    ystore_->UpdateSteps(*s(i));


  return;
}

/*----------------------------------------------------------------------*/
/* compute new direction*/
void STR::INVANA::OptimizerLBFGS::ComputeDirection()
{
  p_->Update(1.0,GetGradientView(),0.0);
  std::vector<double> alpha;

  // loop steps
  for (int i=0; i>-actsize_; i-=numvecs_)
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;

    int ind=0;
    for (int j=numvecs_; j>0; j--)
    {
      MVDotProduct(*(*ystore_)(i-j+1),*(*sstore_)(i-j+1),SolLayoutMapUnique(),&a);
      MVDotProduct(*(*sstore_)(i-j+1),*(*p_)(ind),SolLayoutMapUnique(),&b);
      ind++;
      aa += a;
      bb += b;
    }
    alpha.push_back(1/aa*bb);

    ind=0;
    for (int j=numvecs_; j>0; j--)
    {
      (*p_)(ind)->Update(-1.0*alpha.back(), *(*ystore_)(i-j+1),1.0 );
      ind++;
    }
  }

  // Some scaling of the initial hessian might come in here but has not been proven to be effective
  // altough they say so

  for (int i=-actsize_+1; i<=0; i+=numvecs_)
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;
    double beta=0.0;

    for (int j=0; j<numvecs_; j++)
    {
      MVDotProduct(*(*ystore_)(i+j),*(*sstore_)(i+j),SolLayoutMapUnique(),&a);
      MVDotProduct(*(*ystore_)(i+j),*(*p_)(j),SolLayoutMapUnique(),&b);
      aa += a;
      bb += b;
    }

    beta=1/aa*bb;
    double alphac=alpha.back();
    alpha.pop_back();

    for (int j=0; j<numvecs_; j++)
      (*p_)(j)->Update(alphac-beta, *(*sstore_)(i+j),1.0 );
  }

  // we do minimization not maximization
  p_->Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* print optimization step */
void STR::INVANA::OptimizerLBFGS::PrintOptStep(double tauopt, int numsteps)
{
  double stepincr;
  MVNorm(*step_,SolLayoutMapUnique(),0,&stepincr);

  if (OptProb()->Comm().MyPID()==0)
  {
    printf("STEP %3d | ", runc_);
    printf("Objfunc: %10.8e | ", GetObjFunctValView());
    printf("Gradnorm2: %10.8e | ", convcritc_);
    printf("ErrorIncr: %10.8e | ", error_incr_);
    printf("dp: %10.8e | ", stepincr);
    printf("dt: %10.8e | LSsteps %2d\n", tauopt, numsteps);
    fflush(stdout);
  }

}

/*----------------------------------------------------------------------*/
/* print line search step */
void STR::INVANA::OptimizerLBFGS::PrintLSStep(double tauopt, int numstep)
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
/* print final results*/
void STR::INVANA::OptimizerLBFGS::Summarize()
{
  //std::cout << "the final vector of parameters: " << std::endl;
  //std::cout << *(Matman()->GetParams()) << std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/* Read restart                                               keh 03/14 */
/*----------------------------------------------------------------------*/
void STR::INVANA::OptimizerLBFGS::ReadRestart(int run)
{
  IO::DiscretizationReader reader(OptProb()->Discret(),RestartFromFile(),run);

  if (not OptProb()->Comm().MyPID())
    std::cout << "Reading invana restart from step " << run << " from file: " << RestartFromFile()->FileName() << std::endl;

  //IO::DiscretizationReader reader(discret_, RestartFromFile(),run);
  if (run != reader.ReadInt("run"))
    dserror("Optimization step on file not equal to given step");

  runc_ = run;

  reader.ReadMultiVector(GetSolution(),"optimization_parameters");

  actsize_ = reader.ReadInt("storage_size");

  if (actsize_>0)
  {
    Teuchos::RCP<Epetra_MultiVector> storage = Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(),actsize_,false));
    reader.ReadMultiVector(storage,"sstore");
    for (int i=actsize_-1; i>=0; i--)
      sstore_->UpdateSteps(*(*storage)(i));

    storage->Scale(0.0);
    reader.ReadMultiVector(storage,"ystore");
    for (int i=actsize_-1; i>=0; i--)
      ystore_->UpdateSteps(*(*storage)(i));
  }
  return;
}


/*----------------------------------------------------------------------*/
/* Write restart                                              keh 03/14 */
/*----------------------------------------------------------------------*/
void STR::INVANA::OptimizerLBFGS::WriteRestart()
{
  if (not OptProb()->Comm().MyPID())
    std::cout << "Writing invana restart for step " << runc_ <<  std::endl;

  Writer()->NewStep(runc_, double(runc_));
  Writer()->WriteInt("run", runc_);

  // write vectors with unique gids only
  Teuchos::RCP<Epetra_MultiVector> uniqueparams = Teuchos::rcp(new Epetra_MultiVector(SolLayoutMapUnique(), numvecs_,false));
  LINALG::Export(GetSolutionView(), *uniqueparams);

  Writer()->WriteVector("optimization_parameters", uniqueparams);

  Writer()->WriteInt("storage_size", actsize_);
  //only write vectors if something is there to be written
  if (actsize_>0)
  {
    //write sstore
    Teuchos::RCP<Epetra_MultiVector> uniquestorage = Teuchos::rcp(new Epetra_MultiVector(SolLayoutMapUnique(),actsize_,false));
    for (int i=0; i<actsize_; i++)
      LINALG::Export(*(*sstore_)(-i),*(*uniquestorage)(i));

    Writer()->WriteVector("sstore", uniquestorage);

    // write ystore
    uniquestorage->Scale(0.0);
    for (int i=0; i<actsize_; i++)
      LINALG::Export(*(*ystore_)(-i),*(*uniquestorage)(i));

    Writer()->WriteVector("ystore", uniquestorage);
  }

  return;
}
