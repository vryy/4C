/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_ost.cpp

\brief One-Step-Theta time-integration scheme

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_ost.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid/fluid_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"



/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntOneStepTheta::TimIntOneStepTheta(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      startalgo_(false),
      external_loadsn_(Teuchos::null),
      external_loadsnp_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  FLD::FluidImplicitTimeInt::Init();

  // check, if starting algorithm is desired
  if (numstasteps_ > 0)
  {
    startalgo_ = true;
    if (numstasteps_ > stepmax_) dserror("more steps for starting algorithm than steps overall");
  }

  SetElementTimeParameter();

  CompleteGeneralInit();

  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*-----------------------------------------------------------------------*/
FLD::TimIntOneStepTheta::~TimIntOneStepTheta() { return; }

/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_ == 0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta (theta = %0.2f)   STEP = %4d/%4d \n",
        time_, maxtime_, dta_, theta_, step_, stepmax_);
  }
  return;
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)
  */

  hist_->Update(1.0, *veln_, dta_ * (1.0 - theta_), *accn_, 0.0);


  return;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetStateTimInt()
{
  discret_->SetState("velaf", velnp_);
  if (params_->get<bool>("ost new"))
  {
    if (alefluid_) discret_->SetState("gridvn", gridvn_);
  }
  return;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::CalculateAcceleration(const Teuchos::RCP<const Epetra_Vector> velnp,
    const Teuchos::RCP<const Epetra_Vector> veln, const Teuchos::RCP<const Epetra_Vector> velnm,
    const Teuchos::RCP<const Epetra_Vector> accn, const Teuchos::RCP<Epetra_Vector> accnp)
{
  /*

  Following formulations are for n+1; acceleration values, however, are
  directly stored in vectors at time n (velocity has not yet been updated).

  One-step-Theta:

   acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)

  */

  const double fact1 = 1.0 / (theta_ * dta_);
  const double fact2 = -1.0 / theta_ + 1.0; /* = -1/Theta + 1 */

  accnp->Update(fact1, *velnp, 0.0);
  accnp->Update(-fact1, *veln, 1.0);
  accnp->Update(fact2, *accn, 1.0);

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetGamma(Teuchos::ParameterList& eleparams)
{
  eleparams.set("gamma", theta_);
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::Sep_Multiply()
{
  Sep_->Multiply(false, *velnp_, *fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::OutputofFilteredVel(
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
void FLD::TimIntOneStepTheta::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_time_parameter);
  eleparams.set<int>("Physical Type", physicaltype_);

  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set general element parameters
  eleparams.set("dt", dta_);
  eleparams.set("theta", theta_);
  eleparams.set("omtheta", 1.0 - theta_);

  // set scheme-specific element parameters and vector values
  eleparams.set("total time", time_);

  // full implicit handling of pressure in OST integration
  eleparams.set<int>("ost cont and press", params_->get<int>("ost cont and press"));
  eleparams.set<bool>("ost new", params_->get<bool>("ost new"));

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

void FLD::TimIntOneStepTheta::SetTheta()
{
  // starting algorithm, sets theta = 1.0 for starting steps.
  if (startalgo_)
  {
    // use backward-Euler-type parameter combination
    if (step_ <= numstasteps_)
    {
      if (myrank_ == 0)
      {
        std::cout << "Starting algorithm for OST active. "
                  << "Performing step " << step_ << " of " << numstasteps_
                  << " Backward Euler starting steps" << std::endl;
      }
      theta_ = 1.0;
    }
    else
    {
      // recall original user wish
      theta_ = params_->get<double>("theta");
      // do not enter starting algorithm section in the future
      startalgo_ = false;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
| apply external forces to the fluid                      ghamm 12/2014 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext)
{
  // initialize external force for t_n
  if (step_ <= numstasteps_)
  {
    external_loadsn_ = Teuchos::rcp(new Epetra_Vector(*(*fext)(0)));
    external_loadsnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    external_loads_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  }

  if (external_loadsn_ == Teuchos::null)
    dserror(
        "Starting algorithm for OST missing: Increase number of start steps"
        "or perform initialization of external forces");

  // set external force for t_{n+1}
  external_loadsnp_->Update(1.0, *fext, 0.0);

  // compute interpolated external force at t_{n+\theta}
  external_loads_->Update(theta_, *external_loadsnp_, (1.0 - theta_), *external_loadsn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | output of external forces for restart                     ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::OutputExternalForces()
{
  // call base class
  FLD::FluidImplicitTimeInt::OutputExternalForces();

  if (external_loadsn_ != Teuchos::null)
  {
    output_->WriteVector("fexternal_n", external_loadsn_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | read restart of external forces                           ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::ReadRestart(int step)
{
  // call base class
  FLD::FluidImplicitTimeInt::ReadRestart(step);

  IO::DiscretizationReader reader(discret_, step);
  // check whether external forces were written
  const int have_fexternal = reader.ReadInt("have_fexternal");
  if (have_fexternal != -1)
  {
    external_loadsn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    external_loadsnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    if (step_ > numstasteps_ && params_->get<double>("theta") != 1.0)
    {
      reader.ReadVector(external_loadsn_, "fexternal_n");
      if (have_fexternal != external_loadsn_->GlobalLength())
        dserror("reading of external loads failed");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | update of external forces                                 ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::TimeUpdateExternalForces()
{
  if (external_loadsn_ != Teuchos::null) external_loadsn_->Update(1.0, *external_loadsnp_, 0.0);

  return;
}

/*----------------------------------------------------------------------*|
 | Set Eleparams for turbulence models                          bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::TreatTurbulenceModels(Teuchos::ParameterList& eleparams)
{
  FLD::FluidImplicitTimeInt::TreatTurbulenceModels(eleparams);
  if (reconstructder_)
  {
    FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, velnp_, "velafgrad", alefluid_);
    if (params_->get<bool>("ost new"))
    {
      FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, veln_, "velngrad", alefluid_);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of velocity part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntOneStepTheta::MethodOrderOfAccuracyVel() const
{
  if (theta_ != 0.5)
    return 1;
  else
    return 2;
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of pressure part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntOneStepTheta::MethodOrderOfAccuracyPres() const
{
  if (theta_ != 0.5)
    return 1;
  else
    return 2;
}

/*----------------------------------------------------------------------*
| Return linear error coefficient of velocity           mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
double FLD::TimIntOneStepTheta::MethodLinErrCoeffVel() const
{
  double fac = 0.0;

  if (MethodOrderOfAccuracy() == 1)
    fac = 0.5 - theta_;
  else if (MethodOrderOfAccuracy() == 2)
    fac = -1.0 / 12.0;  // = 1.0/6.0 - 0.5*theta_ with theta_ = 0.5
  else
    dserror("Unknown Order of Accuracy for One Step Theta time integration.");

  return fac;
}
