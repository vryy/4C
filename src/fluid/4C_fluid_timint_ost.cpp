/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for fluids


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_ost.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntOneStepTheta::TimIntOneStepTheta(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      startalgo_(false),
      external_loadsn_(Teuchos::null),
      external_loadsnp_(Teuchos::null)
{
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
    if (numstasteps_ > stepmax_)
      FOUR_C_THROW("more steps for starting algorithm than steps overall");
  }

  set_element_time_parameter();

  CompleteGeneralInit();
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::print_time_step_info()
{
  if (myrank_ == 0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta (theta = %0.2f)   STEP = %4d/%4d \n",
        time_, maxtime_, dta_, theta_, step_, stepmax_);
  }
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::set_old_part_of_righthandside()
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)
  */

  hist_->Update(1.0, *veln_, dta_ * (1.0 - theta_), *accn_, 0.0);
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetStateTimInt()
{
  discret_->set_state("velaf", velnp_);
  if (params_->get<bool>("ost new"))
  {
    if (alefluid_) discret_->set_state("gridvn", gridvn_);
  }
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::calculate_acceleration(const Teuchos::RCP<const Epetra_Vector> velnp,
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
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::SetGamma(Teuchos::ParameterList& eleparams)
{
  eleparams.set("gamma", theta_);
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::Sep_Multiply() { Sep_->Multiply(false, *velnp_, *fsvelaf_); }

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::OutputofFilteredVel(
    Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map();
  Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));

  // get fine scale velocity
  if (scale_sep_ == Inpar::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false, *velnp_, *row_finescaleveltmp);
  else
    FOUR_C_THROW("Unknown separation type!");

  // get filtered or coarse scale velocity
  outvec->Update(1.0, *velnp_, -1.0, *row_finescaleveltmp, 0.0);

  fsoutvec->Update(1.0, *row_finescaleveltmp, 0.0);
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntOneStepTheta::set_element_time_parameter()
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
                  << " Backward Euler starting steps" << '\n';
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
    external_loadsnp_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
    external_loads_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
  }

  if (external_loadsn_ == Teuchos::null)
  {
    FOUR_C_THROW(
        "Starting algorithm for OST missing: Increase number of start steps"
        "or perform initialization of external forces");
  }

  // set external force for t_{n+1}
  external_loadsnp_->Update(1.0, *fext, 0.0);

  // compute interpolated external force at t_{n+\theta}
  external_loads_->Update(theta_, *external_loadsnp_, (1.0 - theta_), *external_loadsn_, 0.0);
}


/*----------------------------------------------------------------------*
 | output of external forces for restart                     ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::output_external_forces()
{
  // call base class
  FLD::FluidImplicitTimeInt::output_external_forces();

  if (external_loadsn_ != Teuchos::null)
  {
    output_->write_vector("fexternal_n", external_loadsn_);
  }
}


/*----------------------------------------------------------------------*
 | read restart of external forces                           ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::read_restart(int step)
{
  // call base class
  FLD::FluidImplicitTimeInt::read_restart(step);

  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::Instance()->InputControlFile(), step);
  // check whether external forces were written
  const int have_fexternal = reader.read_int("have_fexternal");
  if (have_fexternal != -1)
  {
    external_loadsn_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
    external_loadsnp_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
    if (step_ > numstasteps_ && params_->get<double>("theta") != 1.0)
    {
      reader.read_vector(external_loadsn_, "fexternal_n");
      if (have_fexternal != external_loadsn_->GlobalLength())
        FOUR_C_THROW("reading of external loads failed");
    }
  }
}


/*----------------------------------------------------------------------*
 | update of external forces                                 ghamm 12/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::time_update_external_forces()
{
  if (external_loadsn_ != Teuchos::null) external_loadsn_->Update(1.0, *external_loadsnp_, 0.0);
}

/*----------------------------------------------------------------------*|
 | Set Eleparams for turbulence models                          bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntOneStepTheta::treat_turbulence_models(Teuchos::ParameterList& eleparams)
{
  FLD::FluidImplicitTimeInt::treat_turbulence_models(eleparams);
  if (reconstructder_)
  {
    FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, velnp_, "velafgrad", alefluid_);
    if (params_->get<bool>("ost new"))
    {
      FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, veln_, "velngrad", alefluid_);
    }
  }
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of velocity part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntOneStepTheta::method_order_of_accuracy_vel() const
{
  if (theta_ != 0.5)
    return 1;
  else
    return 2;
}

/*----------------------------------------------------------------------*
| Give local order of accuracy of pressure part         mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
int FLD::TimIntOneStepTheta::method_order_of_accuracy_pres() const
{
  if (theta_ != 0.5)
    return 1;
  else
    return 2;
}

/*----------------------------------------------------------------------*
| Return linear error coefficient of velocity           mayr.mt 04/2015 |
*-----------------------------------------------------------------------*/
double FLD::TimIntOneStepTheta::method_lin_err_coeff_vel() const
{
  double fac = 0.0;

  if (method_order_of_accuracy() == 1)
    fac = 0.5 - theta_;
  else if (method_order_of_accuracy() == 2)
    fac = -1.0 / 12.0;  // = 1.0/6.0 - 0.5*theta_ with theta_ = 0.5
  else
    FOUR_C_THROW("Unknown Order of Accuracy for One Step Theta time integration.");

  return fac;
}

FOUR_C_NAMESPACE_CLOSE
