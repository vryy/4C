/*----------------------------------------------------------------------*/
/*! \file
\brief LBFGS optimization algorithm

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#include "optimizer_lbfgs.H"
#include "invana_base.H"
#include "matpar_manager.H"
#include "matpar_manager_elementwise.H"
#include "DcsMatrix.H"
#include "objective_funct.H"

#include "../linalg/linalg_utils.H"
#include "invana_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------*/
/* constructor */
INVANA::OptimizerLBFGS::OptimizerLBFGS(const Teuchos::ParameterList& invp)
    : OptimizerBase(invp),
      initscal_(DRT::INPUT::IntegralValue<bool>(invp, "LBFGSINITSCAL")),
      ssize_(invp.get<int>("SIZESTORAGE")),
      actsize_(0),
      sstore_(Teuchos::null),
      ystore_(Teuchos::null),
      p_(Teuchos::null),
      step_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/* setup algorithm specific stuff */
void INVANA::OptimizerLBFGS::Setup()
{
  actsize_ = 0;

  sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, &SolLayoutMap(), true));
  ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, &SolLayoutMap(), true));

  p_ = Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), 1, true));
  step_ = Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), 1, true));
}

/*----------------------------------------------------------------------*/
/* do the optimization loop*/
void INVANA::OptimizerLBFGS::Integrate()
{
  if (!IsInit()) dserror("OpimizerBase is not inititialzed. Call Init() first");

  int success = 0;

  // solve initially to get quantities:
  Evaluate(GetObjFunctVal(), GetGradient());

  // check gradient by fd:
#if 0
    std::cout << "gradient: " << std::endl;
    PrintDataToScreen(*GetGradient());
    EvaluateGradientFD();

    std::cout << "gradient approx: " << std::endl;
    PrintDataToScreen(*GetGradient());
    exit(0);
#endif

  UpdateGradient();

  // get search direction from gradient:
  p_->Update(-1.0, GetGradientOldView(), 0.0);
  if (DRT::Problem::Instance()->Restart()) ComputeDirection();

  UpdateObjFunctValue();
  error_incr_ = GetObjFunctValView();

  // MVNorm(GetGradientOldView(),SolLayoutMapUnique(),2,&convcritc_);
  GetGradientOldView().Norm2(&convcritc_);

  PrintOptStep(0, 0);

  while (convcritc_ > convtol_ && runc_ < maxiter_)
  {
    double tauopt = 1.0;
    int numsteps = 0;

    // do the line search
    success = EvaluateArmijoRule(&tauopt, &numsteps);
    // -> UpdateSolution() was already called


    if (success == 1)
    {
      std::cout << " Line Search Break Down with current stepsize " << tauopt << std::endl;
      break;
    }

    // get the L2-norm:
    MVNorm(GetGradientView(), SolLayoutMap(), 2, &convcritc_);

    // compute new direction
    StoreVectors();
    ComputeDirection();

    // bring quantities to the next run
    UpdateGradient();
    error_incr_ = GetObjFunctValOldView() - GetObjFunctValView();
    UpdateObjFunctValue();
    runc_++;

    // do some on screen printing
    PrintOptStep(tauopt, numsteps);

    if (restartevry_ and (runc_ % restartevry_ == 0) and runc_) WriteRestart();
  }

  // write final state if not yet written
  if ((restartevry_ == 0) or (runc_ % restartevry_ != 0))
  {
    WriteRestart();
  }

  // append solution to the last output (only if a step was made)
  if (runc_)
  {
    WriteOutput();
  }

  Evaluate(NULL, Teuchos::null);  // set optimization parameters to the optproblem
  Summarize();

  return;
}

/*----------------------------------------------------------------------*/
/* do a line search based on armijo rule */
int INVANA::OptimizerLBFGS::EvaluateArmijoRule(double* tauopt, int* numsteps)
{
  int i = 0;
  int imax = 20;
  double c1 = 1.0e-4;
  // double c2=0.4;
  double tau_max = 1.0e10;
  double tau_min = 1.0e-10;
  double gnorm = 0.0;

  // "last"/"intermediate" values for cubic model
  // these are actually safeguardly set after the first call to the quadratic model
  double tau_l = 0.0;
  double e_l = 0.0;

  int success = 0;

  // safeguard multiplicators
  double blow = 0.1;
  double bhigh = 0.5;

  MVNorm(GetGradientOldView(), SolLayoutMap(), 2, &gnorm);

  double tau_n = std::min(stepsize_, 100.0 / (1 + gnorm));

  // step based on current stepsize
  step_->Update(tau_n, *p_, 0.0);

  while (i < imax && tau_n < tau_max && tau_n > tau_min)
  {
    // make a step
    UpdateSolution(*step_);

    Evaluate(GetObjFunctVal(), GetGradient());

    // check sufficient decrease:
    double dfp_o = 0.0;
    MVDotProduct(GetGradientOldView(), *p_, SolLayoutMap(), &dfp_o);

    if ((GetObjFunctValView() - GetObjFunctValOldView()) < c1 * tau_n * dfp_o)
    {
      double dfp_n = 0.0;
      MVDotProduct(GetGradientView(), *p_, SolLayoutMap(), &dfp_n);
      // check curvature condition:
      if (1)  //(abs(dfp_n)<c2*abs(dfp_o))
      {
        *tauopt = tau_n;
        *numsteps = i + 1;
        return 0;
      }
    }

    // do stepsize prediction based on polynomial models
    if (i == 0)
      success =
          polymod(GetObjFunctValOldView(), dfp_o, tau_n, GetObjFunctValView(), blow, bhigh, tauopt);
    else
      success = polymod(GetObjFunctValOldView(), dfp_o, tau_n, GetObjFunctValView(), blow, bhigh,
          tau_l, e_l, tauopt);

    // repeat if cubic model fails
    if (success == 1)
      success =
          polymod(GetObjFunctValOldView(), dfp_o, tau_n, GetObjFunctValView(), blow, bhigh, tauopt);

    e_l = GetObjFunctValView();
    tau_l = tau_n;
    tau_n = *tauopt;

    // next step based on proposed tauopt;
    step_->Update(tau_n, *p_, 0.0);

    PrintLSStep(tau_l, i);

    UndoUpdateSolution();
    i++;
  }

  return 1;
}

/*----------------------------------------------------------------------*/
/* quadratic model */
int INVANA::OptimizerLBFGS::polymod(
    double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
{
  double lleft = tau_n * blow;
  double lright = tau_n * bhigh;

  *tauopt = -(dfp * tau_n * tau_n) / (2 * (e_n - e_o - dfp * tau_n));
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
/* cubic model model */
int INVANA::OptimizerLBFGS::polymod(double e_o, double dfp, double tau_n, double e_n, double blow,
    double bhigh, double tau_l, double e_l, double* tauopt)
{
  double lleft = tau_n * blow;
  double lright = tau_n * bhigh;

  double a1 = tau_n * tau_n;
  double a2 = tau_n * tau_n * tau_n;
  double a3 = tau_l * tau_l;
  double a4 = tau_l * tau_l * tau_l;

  double deta = a1 * a4 - a2 * a3;

  if (deta < 1.0e-14) return 1;

  double b1 = e_n - (e_o + dfp * tau_n);
  double b2 = e_l - (e_o + dfp * tau_l);

  double c1 = 1 / deta * (a4 * b1 - a2 * b2);
  double c2 = 1 / deta * (-a3 * b1 + a1 * b2);

  double arg = c1 * c1 - 3 * c2 * dfp;

  if (arg < 0.0) return 1;

  *tauopt = (-c1 + sqrt(arg)) / (3 * c2);
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerLBFGS::StoreVectors()
{
  // current storage size
  int actsize = sstore_->NumSteps();
  std::pair<int, int> steps = sstore_->GetSteps();
  int past = steps.first;
  int future = steps.second;

  // increase storage volume if not maximal yet
  if (actsize < ssize_ and
      runc_ != 0)  // in the first call to this we have already a storage from initialization
  {
    sstore_->Resize(past - 1, future, &SolLayoutMap(), false);
    ystore_->Resize(past - 1, future, &SolLayoutMap(), false);
  }

  Epetra_Vector s(SolLayoutMap(), true);

  // push back s
  s.Update(1.0, GetSolutionView(), -1.0, GetSolutionOldView(), 0.0);
  sstore_->UpdateSteps(s);

  // push back y
  s.Update(1.0, GetGradientView(), -1.0, GetGradientOldView(), 0.0);
  ystore_->UpdateSteps(s);


  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerLBFGS::ComputeDirection()
{
  p_->Update(1.0, GetGradientView(), 0.0);

  std::vector<double> alpha;

  // current storage boundary indices
  std::pair<int, int> steps = sstore_->GetSteps();
  int past = steps.first;
  int future = steps.second;

  // loop steps
  for (int i = future; i >= past; i--)
  {
    double a = 0.0;
    double b = 0.0;
    (*sstore_)(i)->Dot((*ystore_)[i], &a);
    p_->Dot((*sstore_)[i], &b);
    alpha.push_back(1 / a * b);

    p_->Update(-1.0 * alpha.back(), (*ystore_)[i], 1.0);
  }

  // initial scaling: see Nocedal, "Numerical Optimization", 2006, p. 178, formula (7.20)
  if (initscal_)
  {
    double nomi = 0.0;
    double denomi = 0.0;
    (*sstore_)(future)->Dot((*ystore_)[future], &nomi);
    (*ystore_)(future)->Dot((*ystore_)[future], &denomi);
    p_->Scale(nomi / denomi);
  }

  for (int i = past; i <= future; i++)
  {
    double a = 0.0;
    double b = 0.0;
    (*sstore_)(i)->Dot((*ystore_)[i], &a);
    p_->Dot((*ystore_)[i], &b);

    double alphac = alpha.back();
    alpha.pop_back();

    p_->Update(alphac - 1 / a * b, (*sstore_)[i], 1.0);
  }

  // we do minimization not maximization
  p_->Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* print optimization step */
void INVANA::OptimizerLBFGS::PrintOptStep(double tauopt, int numsteps)
{
  double stepincr;
  MVNorm(*step_, SolLayoutMap(), 0, &stepincr);

  if (OptProb()->Comm().MyPID() == 0)
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
void INVANA::OptimizerLBFGS::PrintLSStep(double tauopt, int numstep)
{
  if (OptProb()->Comm().MyPID() == 0)
  {
    printf("   LINE SEARCH STEP %3d | ", numstep);
    printf("Objective function: %10.8e | ", GetObjFunctValView());
    printf("stepsize : %10.8e\n", tauopt);
    fflush(stdout);
  }
}

/*----------------------------------------------------------------------*/
/* print final results*/
void INVANA::OptimizerLBFGS::Summarize()
{
  // std::cout << "the final vector of parameters: " << std::endl;
  // std::cout << *(Matman()->GetParams()) << std::endl;
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerLBFGS::ReadRestart(int run)
{
  IO::DiscretizationReader reader(OptProb()->Discret(), RestartFromFile(), run);

  if (not OptProb()->Comm().MyPID())
    std::cout << "Reading invana restart from step " << run
              << " from file: " << RestartFromFile()->FileName() << std::endl;

  runc_ = run;

  reader.ReadMultiVector(GetSolution(), "solution");

  int actsize = reader.ReadInt("storage_size");

  // initialize storage
  sstore_->Resize(-actsize + 1, 0, &SolLayoutMap(), false);
  ystore_->Resize(-actsize + 1, 0, &SolLayoutMap(), false);

  Teuchos::RCP<Epetra_MultiVector> storage =
      Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), actsize, false));

  reader.ReadMultiVector(storage, "sstore");
  for (int i = 0; i < actsize; i++) sstore_->UpdateSteps(*(*storage)(i));

  storage->Scale(0.0);
  reader.ReadMultiVector(storage, "ystore");
  for (int i = 0; i < actsize; i++) ystore_->UpdateSteps(*(*storage)(i));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerLBFGS::WriteRestart()
{
  if (not OptProb()->Comm().MyPID())
    std::cout << "Writing OptimizerLBFGS restart for step " << runc_ << std::endl;

  // initialize a new step
  Writer()->WriteNewStep(runc_, double(runc_));

  // write current solution
  Writer()->WriteNamedVector("solution", GetSolution());

  // prepare for writing storage
  int actsize = sstore_->NumSteps();
  Writer()->WriteNamedInt("storage_size", actsize);

  Teuchos::RCP<Epetra_MultiVector> data =
      Teuchos::rcp(new Epetra_MultiVector(SolLayoutMap(), actsize, false));

  std::pair<int, int> steps = sstore_->GetSteps();
  int past = steps.first;

  // process sstore_
  for (int i = 0; i < actsize; i++) (*data)(i)->Scale(1.0, (*sstore_)[past + i]);

  Writer()->WriteNamedVector("sstore", data);

  // process ystore_
  for (int i = 0; i < actsize; i++) (*data)(i)->Scale(1.0, (*ystore_)[past + i]);

  Writer()->WriteNamedVector("ystore", data);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerLBFGS::WriteOutput()
{
  // map estimate already projected on the elementwise layout
  // Teuchos::RCP<Epetra_MultiVector> result = OptProb()->Matman()->GetRawParams(); //old version
  // reform parameters for final output
  // and do not use metaparametrization                        abirzle 12/2017
  Teuchos::RCP<Epetra_MultiVector> result = OptProb()->Matman()->GetMatParams();

  Writer()->WriteNamedVectors("mean_vb", result);

  // covariance estimation in optimization parameter layout
  double objfuncscal = OptProb()->ObjectiveFunct()->GetScaleFac();
  DcsMatrix covmatrix(sstore_, ystore_, initscal_, true, objfuncscal);

  // get variances in elementwise layout
  Teuchos::RCP<Epetra_MultiVector> stdev = OptProb()->Matman()->GetMatrixDiagonal(covmatrix);

  // make it standard deviations
  double** pointers = stdev->Pointers();
  for (int i = 0; i < stdev->NumVectors(); i++)
  {
    double* const to = pointers[i];
    for (int j = 0; j < stdev->MyLength(); j++)
    {
      // saftey check
      if (to[j] > 0)
        to[j] = sqrt(to[j]);
      else
        std::cout << "WARNING: Your standard deviation is negative. No square calculated."
                  << std::endl;
    }
  }

  Writer()->WriteNamedVectors("stdev_vb", stdev);

  return;
}
