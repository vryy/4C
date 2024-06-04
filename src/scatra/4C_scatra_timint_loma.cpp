/*----------------------------------------------------------------------*/
/*! \file
\brief scatra time integration for loma
\level 2
 *------------------------------------------------------------------------------------------------*/
#include "4C_scatra_timint_loma.hpp"

#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 12/13 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntLoma::ScaTraTimIntLoma(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<CORE::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
      lomaparams_(params),
      initialmass_(0.0),
      thermpressn_(0.0),
      thermpressnp_(0.0),
      thermpressdtn_(0.0),
      thermpressdtnp_(0.0)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                     rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::Init()
{
  // safety check
  if (CORE::UTILS::IntegralValue<int>(*lomaparams_, "SGS_MATERIAL_UPDATE"))
    FOUR_C_THROW(
        "Material update using subgrid-scale temperature currently not supported for loMa "
        "problems. Read remark in file 'scatra_ele_calc_loma.H'!");

  return;
}


/*----------------------------------------------------------------------*
 | setup algorithm                                          rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::Setup()
{
  SetupSplitter();
  return;
}

/*----------------------------------------------------------------------*
 | setup splitter                                          deanda 11/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::SetupSplitter()
{
  // set up a species-temperature splitter (if more than one scalar)
  if (NumScal() > 1)
  {
    splitter_ = Teuchos::rcp(new CORE::LINALG::MapExtractor);
    CORE::LINALG::CreateMapExtractorFromDiscretization(*discret_, NumScal() - 1, *splitter_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | set initial thermodynamic pressure                          vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::set_initial_therm_pressure()
{
  // get thermodynamic pressure from material parameters
  int id = problem_->Materials()->FirstIdByType(CORE::Materials::m_sutherland);
  if (id != -1)  // i.e., Sutherland material found
  {
    const CORE::MAT::PAR::Parameter* mat = problem_->Materials()->ParameterById(id);
    const MAT::PAR::Sutherland* actmat = static_cast<const MAT::PAR::Sutherland*>(mat);

    thermpressn_ = actmat->thermpress_;
  }
  else
  {
    FOUR_C_THROW(
        "No Sutherland material found for initial setting of "
        "thermodynamic pressure!");
  }

  // initialize also value at n+1
  // (computed if not constant, otherwise prescribed value remaining)
  thermpressnp_ = thermpressn_;

  // initialize time derivative of thermodynamic pressure at n+1 and n
  // (computed if not constant, otherwise remaining zero)
  thermpressdtnp_ = 0.0;
  thermpressdtn_ = 0.0;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  // -> For constant thermodynamic pressure, this is done here once and
  // for all simulation time.
  compute_therm_pressure_intermediate_values();

  return;
}  // SCATRA::ScaTraTimIntLoma::set_initial_therm_pressure


/*----------------------------------------------------------------------*
 | compute initial total mass in domain                        vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::ComputeInitialMass()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->set_state("phinp", phin_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_total_and_mean_scalars, eleparams);
  // inverted scalar values are required here
  eleparams.set("inverting", true);
  eleparams.set("calc_grad_phi", false);

  // evaluate integral of inverse temperature
  Teuchos::RCP<CORE::LINALG::SerialDenseVector> scalars =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // compute initial mass times gas constant: R*M_0 = int(1/T_0)*tp
  initialmass_ = (*scalars)[0] * thermpressn_;

  // print out initial total mass
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
    std::cout << "Initial total mass in domain (times gas constant): " << initialmass_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
  }

  return;
}  // SCATRA::ScaTraTimIntLoma::ComputeInitialMass


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure from mass conservation       vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::compute_therm_pressure_from_mass_cons()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->set_state("phinp", phinp_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_total_and_mean_scalars, eleparams);
  // inverted scalar values are required here
  eleparams.set("inverting", true);
  eleparams.set("calc_grad_phi", false);

  // evaluate integral of inverse temperature
  Teuchos::RCP<CORE::LINALG::SerialDenseVector> scalars =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // compute thermodynamic pressure: tp = R*M_0/int(1/T)
  thermpressnp_ = initialmass_ / (*scalars)[0];

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
    std::cout << "Thermodynamic pressure from mass conservation: " << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  compute_therm_pressure_time_derivative();

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  compute_therm_pressure_intermediate_values();

  return;
}  // SCATRA::ScaTraTimIntLoma::compute_therm_pressure_from_mass_cons


/*----------------------------------------------------------------------*
 | add parameters depending on the problem              rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntLoma::add_problem_specific_parameters_and_vectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  add_therm_press_to_parameter_list(params);
  return;
}

FOUR_C_NAMESPACE_CLOSE
