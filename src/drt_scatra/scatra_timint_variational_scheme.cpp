/*----------------------------------------------------------------------*/
/*!

\brief  Time integration for variational formulation problems using different schemes

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_variational_scheme.H"
// To read initial conditions from input file
#include "../drt_lib/drt_globalproblem.H"
// To read parameters from input file
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                   deanda 10/17 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntVariationalOST::TimIntVariationalOST(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output, probnum),
      TimIntVariational(actdis, solver, params, extraparams, output, probnum),
      TimIntOneStepTheta(actdis, solver, params, extraparams, output, probnum),
      semImplicitFunctional_(
          DRT::INPUT::IntegralValue<int>(params->sublist("VARIATIONAL"), "SEMIMPLICITFUNCTIONAL"))
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                            deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntVariational::Init();
  TimIntOneStepTheta::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                            deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  // back-up of the initial condition
  TimIntVariational::Setup();
  TimIntOneStepTheta::Setup();
  return;
}

void SCATRA::TimIntVariationalOST::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  // call base class routine
  TimIntVariational::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
  discret_->SetState("phin", phin_);
  TimIntOneStepTheta::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);

  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                 deanda 10/17 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntVariationalOST::~TimIntVariationalOST() { return; }


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                         deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::Update(const int num)
{
  TimIntVariational::Update(num);
  TimIntOneStepTheta::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | explicit predictor for nonlinear solver                 deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::ExplicitPredictor() const
{
  // call base class routine
  TimIntVariational::ExplicitPredictor();
  TimIntOneStepTheta::ExplicitPredictor();

  return;
}

/*----------------------------------------------------------------------*
 | explicit predictor for nonlinear solver                 deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::PrepareFirstTimeStep()
{
  TimIntVariational::ApplyDirichletBC(time_, phin_, Teuchos::null);
  TimIntVariational::PrepareFirstTimeStep();

  return;
}  // SCATRA::TimIntVariationalOST::PrepareFirstTimeStep

/*----------------------------------------------------------------------*
 | add parameters depending on the problem                  deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariationalOST::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  TimIntVariational::AddProblemSpecificParametersAndVectors(params);

  // Defines the type of scheme to use in the variational functional: Implicit or semi-implicit
  params.set<bool>("Is_semImplicit_Functional", semImplicitFunctional_);
  return;
}
