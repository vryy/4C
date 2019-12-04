/*-----------------------------------------------------------*/
/*! \file

\brief BDF-2 time-integration scheme

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_bdf2.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntBDF2::TimIntBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid), theta_(1.0)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  FLD::FluidImplicitTimeInt::Init();

  // check, if starting algorithm is desired
  if (numstasteps_ > 0)
    dserror("no starting algorithm supported for schemes other than af-gen-alpha");

  SetElementTimeParameter();

  CompleteGeneralInit();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntBDF2::~TimIntBDF2() { return; }

/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::PrintTimeStepInfo()
{
  if (myrank_ == 0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n", time_,
        maxtime_, dta_, step_, stepmax_);
  }
  return;
}

/*----------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetTheta()
{
  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt

  if (step_ > 1)
    theta_ = (dta_ + dtp_) / (2.0 * dta_ + dtp_);
  else
  {
    // use backward Euler for the first time step
    velnm_->Update(1.0, *veln_, 0.0);  // results in hist_ = veln_
    theta_ = 1.0;
  }

  return;
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetOldPartOfRighthandside()
{
  /*
     BDF2: for constant time step:

                   mom: hist_ = 4/3 veln_  - 1/3 velnm_
                  (con: hist_ = 4/3 densn_ - 1/3 densnm_)

  */

  hist_->Update(4. / 3., *veln_, -1. / 3., *velnm_, 0.0);

  return;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetStateTimInt()
{
  discret_->SetState("velaf", velnp_);

  return;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::CalculateAcceleration(const Teuchos::RCP<const Epetra_Vector> velnp,
    const Teuchos::RCP<const Epetra_Vector> veln, const Teuchos::RCP<const Epetra_Vector> velnm,
    const Teuchos::RCP<const Epetra_Vector> accn, const Teuchos::RCP<Epetra_Vector> accnp)
{
  /*

  BDF2:

                 2*dt(n)+dt(n-1)                  dt(n)+dt(n-1)
   acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
               dt(n)*[dt(n)+dt(n-1)]              dt(n)*dt(n-1)

                       dt(n)
             + ----------------------- vel(n-1)
               dt(n-1)*[dt(n)+dt(n-1)]

  */

  if (dta_ * dtp_ < EPS15) dserror("Zero time step size!!!!!");
  const double sum = dta_ + dtp_;

  accnp->Update((2.0 * dta_ + dtp_) / (dta_ * sum), *velnp, -sum / (dta_ * dtp_), *veln, 0.0);
  accnp->Update(dta_ / (dtp_ * sum), *velnm, 1.0);

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetGamma(Teuchos::ParameterList& eleparams)
{
  eleparams.set("gamma", 1.0);
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::Sep_Multiply()
{
  Sep_->Multiply(false, *velnp_, *fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntBDF2::OutputofFilteredVel(
    Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));

  // get fine scale velocity
  if (scale_sep_ == INPAR::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false, *velnp_, *row_finescaleveltmp);
  else
    dserror("Unknown separation type!");

  // get filtered or coarse scale velocity
  outvec->Update(1.0, *velnp_, -1.0, *row_finescaleveltmp, 0.0);

  fsoutvec->Update(1.0, *row_finescaleveltmp, 0.0);

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntBDF2::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_time_parameter);

  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set general element parameters
  eleparams.set("dt", dta_);
  eleparams.set("theta", theta_);
  eleparams.set("omtheta", 0.0);

  // set scheme-specific element parameters and vector values
  eleparams.set("total time", time_);


  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| Return linear error coefficient of velocity             mayr.mt 12/13 |
*-----------------------------------------------------------------------*/
double FLD::TimIntBDF2::MethodLinErrCoeffVel() const
{
  double nominator = (dta_ + dtp_) * (dta_ + dtp_);
  double denominator = 6 * dta_ * (2 * dta_ + dtp_);

  return nominator / denominator;
}
