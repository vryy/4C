/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_ptc.cpp

\brief pseudo transient solution method

\maintainer Michael Hiermeier

\date Oct 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_ptc.H"   // class definition
#include "str_nln_solver_utils.H"
#include "str_timint_base.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::PseudoTransient::PseudoTransient()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::PseudoTransient::Setup()
{
  CheckInit();

  // setup the nox parameter list for a pseudo transient solution method
  SetPseudoTransientParams();

  // Call the Setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  Nox::Setup();

  if (not IsSetup())
    dserror("issetup_ should be \"true\" at this point!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::PseudoTransient::SetPseudoTransientParams()
{
  CheckInit();

  /* get the nox parameter list and set the necessary parameters for a
   * pseudo transient solution procedure */
  Teuchos::ParameterList& pnox = DataSDyn().GetMutableNoxParams();

  // ---------------------------------------------------------------------------
  // Set-up the pseudo transient method
  // ---------------------------------------------------------------------------
  // Non-linear solver
  pnox.set("Nonlinear Solver","Pseudo Transient");

  // Direction
  Teuchos::ParameterList& pdir = pnox.sublist("Direction",true);
  pdir.set("Method","Newton");

  // Line Search
  Teuchos::ParameterList& plinesearch = pnox.sublist("Line Search",true);

  // Line Search/Full Step (line search is deactivated)
  Teuchos::ParameterList& pfullstep = plinesearch.sublist("Full Step",true);
  // check if the default value is set
  double fullstep = pfullstep.get<double>("Full Step");
  if (fullstep!=1.0)
  {
    std::string markerline;
    markerline.assign(40,'!');
    std::cout << markerline << std::endl
        << "WARNING: The Full Step length is " << fullstep
        << " (default=1.0)"<< std::endl
        << markerline << std::endl;
  }
  /* The following parameters create a NOX::NLN::Solver::PseudoTransient
   * solver which is equivalent to the old BACI implementation.
   *
   * If you are keen on using the new features, please use the corresponding
   * input section "STUCT NOX/Pseudo Transient" in your input file. */
  Teuchos::ParameterList& pptc = pnox.sublist("Pseudo Transient");

  pptc.set<double>("deltaInit",DataSDyn().GetInitialPTCPseudoTimeStep());
  pptc.set<double>("deltaMax",std::numeric_limits<double>::max());
  pptc.set<double>("deltaMin",0.0);
  pptc.set<int>("Maximum Number of Pseudo-Transient Iterations",
      (DataSDyn().GetIterMax()+1));
  pptc.set<std::string>("Time Step Control","SER");
  pptc.set<std::string>("Norm Type for TSC","Max Norm");
  pptc.set<std::string>("Scaling Type","Identity");

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the dat-file instead.
   */
  if (not IsXMLStatusTestFile(DataSDyn().GetNoxParams().sublist("Status Test")))
  {
    std::set<enum NOX::NLN::StatusTest::QuantityType> qtypes;
    CreateQuantityTypes(qtypes,DataSDyn());
    SetStatusTestParams(DataSDyn().GetMutableNoxParams().sublist("Status Test"),
        DataSDyn(),qtypes);
  }

  return;
}
