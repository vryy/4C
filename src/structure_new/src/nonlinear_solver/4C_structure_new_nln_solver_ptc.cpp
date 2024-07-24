/*-----------------------------------------------------------*/
/*! \file

\brief pseudo transient solution method


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_solver_ptc.hpp"  // class definition

#include "4C_structure_new_nln_solver_utils.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Nln::SOLVER::PseudoTransient::PseudoTransient()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Nln::SOLVER::PseudoTransient::setup()
{
  check_init();

  // setup the nox parameter list for a pseudo transient solution method
  set_pseudo_transient_params();

  // Call the setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  Nox::setup();

  FOUR_C_ASSERT(is_setup(), "issetup_ should be \"true\" at this point!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Nln::SOLVER::PseudoTransient::set_pseudo_transient_params()
{
  check_init();

  /* get the nox parameter list and set the necessary parameters for a
   * pseudo transient solution procedure */
  Teuchos::ParameterList& pnox = data_sdyn().get_nox_params();

  // ---------------------------------------------------------------------------
  // Set-up the pseudo transient method
  // ---------------------------------------------------------------------------
  // Non-linear solver
  pnox.set("Nonlinear Solver", "Pseudo Transient");

  // Direction
  Teuchos::ParameterList& pdir = pnox.sublist("Direction", true);
  pdir.set("Method", "Newton");

  // Line Search
  Teuchos::ParameterList& plinesearch = pnox.sublist("Line Search", true);

  // Line Search/Full Step (line search is deactivated)
  Teuchos::ParameterList& pfullstep = plinesearch.sublist("Full Step", true);
  // check if the default value is set
  double fullstep = pfullstep.get<double>("Full Step");
  if (fullstep != 1.0)
  {
    std::string markerline;
    markerline.assign(40, '!');
    std::cout << markerline << std::endl
              << "WARNING: The Full Step length is " << fullstep << " (default=1.0)" << std::endl
              << markerline << std::endl;
  }
  /* The following parameters create a NOX::Nln::Solver::PseudoTransient
   * solver which is equivalent to the old 4C implementation.
   *
   * If you are keen on using the new features, please use the corresponding
   * input section "STRUCT NOX/Pseudo Transient" in your input file. */
  Teuchos::ParameterList& pptc = pnox.sublist("Pseudo Transient");

  pptc.set<double>("deltaInit", data_sdyn().get_initial_ptc_pseudo_time_step());
  pptc.set<double>("deltaMax", std::numeric_limits<double>::max());
  pptc.set<double>("deltaMin", 0.0);
  pptc.set<int>("Maximum Number of Pseudo-Transient Iterations", (data_sdyn().get_iter_max() + 1));
  pptc.set<std::string>("Time Step Control", "SER");
  pptc.set<double>("SER_alpha", 1.0);
  pptc.set<double>("ScalingFactor", 1.0);
  pptc.set<std::string>("Norm Type for TSC", "Max Norm");
  pptc.set<std::string>("Build scaling operator", "every timestep");
  pptc.set<std::string>("Scaling Type", "Identity");

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the dat-file instead.
   */
  if (not is_xml_status_test_file(data_sdyn().get_nox_params().sublist("Status Test")))
  {
    std::set<enum NOX::Nln::StatusTest::QuantityType> qtypes;
    create_quantity_types(qtypes, data_sdyn());
    set_status_test_params(
        data_sdyn().get_nox_params().sublist("Status Test"), data_sdyn(), qtypes);
  }
}

FOUR_C_NAMESPACE_CLOSE
