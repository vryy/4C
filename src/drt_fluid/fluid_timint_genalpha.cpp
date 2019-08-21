/*-----------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time-integration scheme

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/


#include "fluid_timint_genalpha.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntGenAlpha::TimIntGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      alphaM_(params_->get<double>("alpha_M")),
      alphaF_(params_->get<double>("alpha_F")),
      gamma_(params_->get<double>("gamma")),
      startalgo_(false)
// af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
// (may be reset below when starting algorithm is used)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  FLD::FluidImplicitTimeInt::Init();

  // starting algorithm only for af-generalized-alpha so far
  // -> check for time-integration scheme and reasonable number of steps
  if (numstasteps_ > 0)
  {
    if (timealgo_ != INPAR::FLUID::timeint_afgenalpha)
      dserror("no starting algorithm supported for schemes other than af-gen-alpha");
    else
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
FLD::TimIntGenAlpha::~TimIntGenAlpha() { return; }


/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::PrintTimeStepInfo()
{
  if (myrank_ == 0)
  {
    switch (timealgo_)
    {
      case INPAR::FLUID::timeint_afgenalpha:
        printf(
            "TIME: %11.4E/%11.4E  DT = %11.4E  Af-Generalized-Alpha (gamma = %0.2f, alphaF = "
            "%0.2f, alphaM = %0.2f) STEP = %4d/%4d \n",
            time_, maxtime_, dta_, gamma_, alphaF_, alphaM_, step_, stepmax_);
        break;
      case INPAR::FLUID::timeint_npgenalpha:
        printf(
            "TIME: %11.4E/%11.4E  DT = %11.4E  Np-Generalized-Alpha (gamma = %0.2f, alphaF = "
            "%0.2f, alphaM = %0.2f) STEP = %4d/%4d \n",
            time_, maxtime_, dta_, gamma_, alphaF_, alphaM_, step_, stepmax_);
        break;
      default:
        dserror("parameter out of range: IOP\n");
        break;
    } /* end of switch(timealgo) */
  }
  return;
}


/*----------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetTheta()
{
  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  // starting algorithm
  if (startalgo_)
  {
    // use backward-Euler-type parameter combination
    if (step_ <= numstasteps_)
    {
      if (myrank_ == 0)
      {
        std::cout << "Starting algorithm for Af_GenAlpha active. "
                  << "Performing step " << step_ << " of " << numstasteps_
                  << " Backward Euler starting steps" << std::endl;
      }
      alphaM_ = 1.0;
      alphaF_ = 1.0;
      gamma_ = 1.0;
    }
    else
    {
      // recall original user wish
      alphaM_ = params_->get<double>("alpha_M");
      alphaF_ = params_->get<double>("alpha_F");
      gamma_ = params_->get<double>("gamma");
      // do not enter starting algorithm section in the future
      startalgo_ = false;
    }
  }

  // compute "pseudo-theta" for af-generalized-alpha scheme
  theta_ = alphaF_ * gamma_ / alphaM_;

  return;
}


/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetOldPartOfRighthandside()
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

      af-generalized-alpha:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)

  */

  hist_->PutScalar(0.0);

  return;
}


/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::GenAlphaUpdateAcceleration()
{
  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // compute factors
  const double fact1 = 1.0 / (gamma_ * dta_);
  const double fact2 = 1.0 - (1.0 / gamma_);

  // consider both velocity and pressure degrees of freedom in case of
  // artificial compressibility or weakly_compressible or
  // extract and update only velocity degrees of freedom, since in
  // low-Mach-number flow, 'pressure' components are used to store
  // temporal derivatives of scalar/temperature values
  if (physicaltype_ == INPAR::FLUID::artcomp or
      physicaltype_ == INPAR::FLUID::weakly_compressible or
      physicaltype_ == INPAR::FLUID::weakly_compressible_stokes)
  {
    accnp_->Update(fact2, *accn_, 0.0);
    accnp_->Update(fact1, *velnp_, -fact1, *veln_, 1.0);
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn = velpressplitter_->ExtractOtherVector(accn_);
    Teuchos::RCP<Epetra_Vector> onlyveln = velpressplitter_->ExtractOtherVector(veln_);
    Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_->ExtractOtherVector(velnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

    onlyaccnp->Update(fact2, *onlyaccn, 0.0);
    onlyaccnp->Update(fact1, *onlyvelnp, -fact1, *onlyveln, 1.0);

    // copy back into global vector
    LINALG::Export(*onlyaccnp, *accnp_);
  }

}  // TimIntGenAlpha::GenAlphaUpdateAcceleration


/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::GenAlphaIntermediateValues()
{
  // set intermediate values for acceleration and potential temporal
  // derivatives
  //
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)

  // consider both velocity and pressure degrees of freedom in case of
  // artificial compressibility or weakly_compressible or
  // extract and update only velocity degrees of freedom, since in
  // low-Mach-number flow, 'pressure' components are used to store
  // temporal derivatives of scalar/temperature values
  if (physicaltype_ == INPAR::FLUID::artcomp or
      physicaltype_ == INPAR::FLUID::weakly_compressible or
      physicaltype_ == INPAR::FLUID::weakly_compressible_stokes)
  {
    accam_->Update((alphaM_), *accnp_, (1.0 - alphaM_), *accn_, 0.0);
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn = velpressplitter_->ExtractOtherVector(accn_);
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_->ExtractOtherVector(accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_), *onlyaccnp, (1.0 - alphaM_), *onlyaccn, 0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam, *accam_);
  }

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  velaf_->Update((alphaF_), *velnp_, (1.0 - alphaF_), *veln_, 0.0);

}  // TimIntGenAlpha::GenAlphaIntermediateValues


/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha             |
 | for given vectors at n and n+1                           ghamm 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::GenAlphaIntermediateValues(
    Teuchos::RCP<Epetra_Vector>& vecnp, Teuchos::RCP<Epetra_Vector>& vecn)
{
  // compute intermediate values for given vectors
  //
  //       n+alphaM
  //    vec         = alpha_M * vecnp     + (1-alpha_M) *  vecn
  //
  //       n+alphaF
  //    vec         = alpha_F * vecnp     + (1-alpha_F) *  vecn

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> vecmap = Teuchos::rcp(new Epetra_Map(vecnp->Map().NumGlobalElements(),
      vecnp->Map().NumMyElements(), vecnp->Map().MyGlobalElements(), 0, vecnp->Map().Comm()));

  Teuchos::RCP<Epetra_Vector> vecam = LINALG::CreateVector(*vecmap, true);
  vecam->Update((alphaM_), *vecnp, (1.0 - alphaM_), *vecn, 0.0);

  Teuchos::RCP<Epetra_Vector> vecaf = LINALG::CreateVector(*vecmap, true);
  vecaf->Update((alphaF_), *vecnp, (1.0 - alphaF_), *vecn, 0.0);

  // store computed intermediate values in given vectors
  vecnp = vecaf;
  vecn = vecam;

  return;
}  // TimIntGenAlpha::GenAlphaIntermediateValues


/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetStateTimInt()
{
  discret_->SetState("velaf", velaf_);
  discret_->SetState("velam", velam_);
  if (timealgo_ == INPAR::FLUID::timeint_npgenalpha) discret_->SetState("velnp", velnp_);

  return;
}

/*----------------------------------------------------------------------*|
 | Set Eleparams for turbulence models                          bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::TreatTurbulenceModels(Teuchos::ParameterList& eleparams)
{
  FLD::FluidImplicitTimeInt::TreatTurbulenceModels(eleparams);
  if (reconstructder_)
    FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, velaf_, "velafgrad", alefluid_);
  return;
}

/*----------------------------------------------------------------------*
| return alphaF_                                               bk 12/13 |
*-----------------------------------------------------------------------*/
double FLD::TimIntGenAlpha::SetTimeFac() { return alphaF_; }


/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::CalculateAcceleration(const Teuchos::RCP<const Epetra_Vector> velnp,
    const Teuchos::RCP<const Epetra_Vector> veln, const Teuchos::RCP<const Epetra_Vector> velnm,
    const Teuchos::RCP<const Epetra_Vector> accn, const Teuchos::RCP<Epetra_Vector> accnp)
{
  // do nothing: new acceleration is calculated at beginning of next time step

  return;
}


/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetGamma(Teuchos::ParameterList& eleparams)
{
  eleparams.set("gamma", gamma_);
  return;
}


/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::Sep_Multiply()
{
  Sep_->Multiply(false, *velaf_, *fsvelaf_);
  return;
}


/*----------------------------------------------------------------------*
| update velaf_                                                bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::UpdateVelafGenAlpha()
{
  velaf_->Update((alphaF_), *velnp_, (1.0 - alphaF_), *veln_, 0.0);
  return;
}


/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::OutputofFilteredVel(
    Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));

  // get fine scale velocity
  if (scale_sep_ == INPAR::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false, *velaf_, *row_finescaleveltmp);
  else
    dserror("Unknown separation type!");

  // get filtered or coarse scale velocity
  outvec->Update(1.0, *velaf_, -1.0, *row_finescaleveltmp, 0.0);

  fsoutvec->Update(1.0, *row_finescaleveltmp, 0.0);

  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntGenAlpha::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_time_parameter);
  eleparams.set<int>("Physical Type", physicaltype_);

  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);
  // set general element parameters
  eleparams.set("dt", dta_);
  eleparams.set("theta", theta_);
  eleparams.set("omtheta", 0.0);

  // set scheme-specific element parameters and vector values
  if (time_ > 0.0)
  {
    eleparams.set("total time", time_ - (1 - alphaF_) * dta_);
  }
  else
  {
    eleparams.set("total time", time_);
  }
  eleparams.set("alphaF", alphaF_);
  eleparams.set("alphaM", alphaM_);
  eleparams.set("gamma", gamma_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| extrapolate end point                                        bk 12/13 |
*-----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::TimIntGenAlpha::ExtrapolateEndPoint(
    Teuchos::RCP<Epetra_Vector> vecn, Teuchos::RCP<Epetra_Vector> vecm)
{
  Teuchos::RCP<Epetra_Vector> vecnp = FluidImplicitTimeInt::ExtrapolateEndPoint(vecn, vecm);

  // For gen-alpha extrapolate mid-point quantities to end-point.
  // Otherwise, equilibrium time level is already end-point.

  vecnp->Update((alphaF_ - 1.0) / alphaF_, *vecn, 1.0 / alphaF_);

  return vecnp;
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of velocity part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntGenAlpha::MethodOrderOfAccuracyVel() const
{
  if (fabs(gamma_ - 0.5 - alphaM_ + alphaF_) < 1.0e-6)
    return 2;
  else
    return 1;
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of pressure part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntGenAlpha::MethodOrderOfAccuracyPres() const
{
  if (fabs(gamma_ - 0.5 - alphaM_ + alphaF_) < 1.0e-6)
    return 2;
  else
    return 1;
}

/*----------------------------------------------------------------------*
| Return linear error coefficient of velocity           mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
double FLD::TimIntGenAlpha::MethodLinErrCoeffVel() const
{
  double fac = 0.0;

  if (MethodOrderOfAccuracy() == 1)
    fac = 0.5 - gamma_;
  else if (MethodOrderOfAccuracy() == 2)
    fac = 1.0 / 6.0 - 0.5 * gamma_;
  else
    dserror("Unknown Order of Accuracy for Gen-Alpha time integration.");

  return fac;
}
