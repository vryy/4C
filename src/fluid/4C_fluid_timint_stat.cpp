/*-----------------------------------------------------------*/
/*! \file

\brief Basic fluid driver for stationary problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_stat.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_turbulence_boxfilter.hpp"
#include "4C_fluid_turbulence_dyn_smag.hpp"
#include "4C_fluid_turbulence_dyn_vreman.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fluid_volumetric_surfaceFlow_condition.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntStationary::TimIntStationary(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntStationary::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  FLD::FluidImplicitTimeInt::Init();

  // check, if starting algorithm is desired
  if (numstasteps_ > 0)
    FOUR_C_THROW("no starting algorithm supported for schemes other than af-gen-alpha");

  set_element_time_parameter();

  CompleteGeneralInit();
  return;
}



/*----------------------------------------------------------------------*
 | use TimeLoop() to start stationary problem                  bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntStationary::TimeLoop()
{
  solve_stationary_problem();
  return;
}


/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationary::set_old_part_of_righthandside()
{
  /*
     Stationary:

                   mom: hist_ = 0.0
                   (con: hist_ = 0.0)
  */

  hist_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 | solve stationary fluid problem                              gjb 10/07|
 *----------------------------------------------------------------------*/
void FLD::TimIntStationary::solve_stationary_problem()
{
  // time measurement: time loop (stationary) --- start TimeMonitor tm2
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  while (NotFinished())
  {
    // -------------------------------------------------------------------
    //              set (pseudo-)time-dependent parameters
    // -------------------------------------------------------------------
    increment_time_and_step();

    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_ == 0)
    {
      printf("Stationary Fluid Solver - STEP = %4d/%4d \n", step_, stepmax_);
    }

    set_element_time_parameter();

    // -------------------------------------------------------------------
    //         evaluate Dirichlet and Neumann boundary conditions
    // -------------------------------------------------------------------
    set_dirichlet_neumann_bc();

    // -------------------------------------------------------------------
    //           preparation of AVM3-based scale separation
    // -------------------------------------------------------------------
    if (step_ == 1 and fssgv_ != INPAR::FLUID::no_fssgv) AVM3Preparation();

    // -------------------------------------------------------------------
    //                     solve equation system
    // -------------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //                        compute flow rates
    // -------------------------------------------------------------------
    ComputeFlowRates();

    // -------------------------------------------------------------------
    // evaluate divergence u
    // -------------------------------------------------------------------
    EvaluateDivU();


    // Evaluate integral error
    evaluate_error_compared_to_analytical_sol();



    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  }  // end of time loop

}  // TimIntStationary::solve_stationary_problem

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationary::SetStateTimInt()
{
  discret_->set_state("velaf", velnp_);
  return;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationary::calculate_acceleration(const Teuchos::RCP<const Epetra_Vector> velnp,
    const Teuchos::RCP<const Epetra_Vector> veln, const Teuchos::RCP<const Epetra_Vector> velnm,
    const Teuchos::RCP<const Epetra_Vector> accn, const Teuchos::RCP<Epetra_Vector> accnp)
{
  accnp->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationary::SetGamma(Teuchos::ParameterList& eleparams)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationary::Sep_Multiply()
{
  Sep_->Multiply(false, *velnp_, *fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntStationary::OutputofFilteredVel(
    Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  // no output since subgrid-scale modeling does not make sense for stationary problems!!!
  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntStationary::set_element_time_parameter()
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
  eleparams.set("total time", time_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| return time integration factor                               bk 12/13 |
*-----------------------------------------------------------------------*/
double FLD::TimIntStationary::TimIntParam() const
{
  double retval = 0.0;
  // no FSI with stationary time integrator
  FOUR_C_THROW("FSI does not allow a stationary time integrator.");
  return retval;
}

/*----------------------------------------------------------------------*|
 | Set Eleparams for turbulence models                          bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntStationary::treat_turbulence_models(Teuchos::ParameterList& eleparams)
{
  FLD::FluidImplicitTimeInt::treat_turbulence_models(eleparams);
  if (reconstructder_)
    FLD::UTILS::ProjectGradientAndSetParam(discret_, eleparams, velnp_, "velafgrad", alefluid_);
  return;
}

FOUR_C_NAMESPACE_CLOSE
