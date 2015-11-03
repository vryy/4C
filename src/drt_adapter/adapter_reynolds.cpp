/*--------------------------------------------------------------------------*/
/*!
\file adapter_reynolds.cpp

\brief Reynolds field base algorithm

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "adapter_reynolds.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_solver.H"

#include "../drt_scatra/scatra_resulttest.H"

#include "../drt_scatra/scatra_timint_stat.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ReynoldsBaseAlgorithm::Setup(
    const Teuchos::ParameterList&   prbdyn,         ///< parameter list for global problem
    const Teuchos::ParameterList&   reynoldsdyn,    ///< parameter list for Reynolds subproblem
    const Teuchos::ParameterList&   solverparams,   ///< parameter list for Reynolds solver
    const std::string&              disname,        ///< name of Reynolds discretization
    const bool                      isale           ///< ALE flag
    )
{
  // setup Reynolds algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter list? -> no default paramter possible any more
  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(solverparams,
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  Teuchos::RCP<Teuchos::ParameterList> reynoldstimeparams = Teuchos::rcp(new Teuchos::ParameterList(reynoldsdyn));

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  // the default time step size
  reynoldstimeparams->set<double>   ("TIMESTEP"    ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  reynoldstimeparams->set<double>   ("MAXTIME"     ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  reynoldstimeparams->set<int>      ("NUMSTEP"     ,prbdyn.get<int>("NUMSTEP"));
  // restart
  reynoldstimeparams->set           ("RESTARTEVRY" ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  reynoldstimeparams->set           ("UPRES"       ,prbdyn.get<int>("UPRES"));

  // -------------------------------------------------------------------
  // list for extra parameters
  // (put here everything that is not available in reynoldsdyn or its sublists)
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> extraparams
    = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale",isale);

  // ------------------------------------get also fluid turbulence sublist
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  extraparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  extraparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

    // create instance of time integration class (call the constructor)
  reynolds_ = Teuchos::rcp(
      new SCATRA::TimIntStationary(actdis, solver, reynoldstimeparams,
          extraparams, output));

  reynolds_->Init();
  // initialize algorithm for specific time-integration scheme

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ReynoldsBaseAlgorithm::CreateReynoldsFieldTest()
{
  return Teuchos::rcp(new SCATRA::ScaTraResultTest(reynolds_));
}
