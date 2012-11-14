/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo.cpp

\brief Thermo field adapter

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                             bborn 08/09 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                 bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "adapter_thermo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_thermo/thrtimint_statics.H"
#include "../drt_thermo/thrtimint_genalpha.H"
#include "../drt_thermo/thrtimint_ost.H"
#include "../drt_thermo/thrtimint_expleuler.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for ThermoBaseAlgorithm:
#include "../drt_inpar/inpar_thermo.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
ADAPTER::Thermo::~Thermo()
{
}

/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
ADAPTER::ThermoBaseAlgorithm::ThermoBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  SetupThermo(prbdyn);
}

/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
ADAPTER::ThermoBaseAlgorithm::~ThermoBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoBaseAlgorithm::SetupThermo(const Teuchos::ParameterList& prbdyn)
{
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();

  // major switch to different time integrators
  INPAR::THR::DynamicType timinttype = DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn,"DYNAMICTYP");
  switch (timinttype)
  {
  case INPAR::THR::dyna_statics :
  case INPAR::THR::dyna_onesteptheta :
  case INPAR::THR::dyna_gemm :
  case INPAR::THR::dyna_genalpha :
  case INPAR::THR::dyna_expleuler :
    SetupTimInt(prbdyn, timinttype);   // <-- here is the show
    break;
  default :
    dserror("unknown time integration scheme '%s'", tdyn.get<std::string>("DYNAMICTYP").c_str());
    break;
  }

}

/*----------------------------------------------------------------------*
 | setup of thermal time integration                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoBaseAlgorithm::SetupTimInt(
  const Teuchos::ParameterList& prbdyn,
  INPAR::THR::DynamicType timinttype
  )
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t
    = Teuchos::TimeMonitor::getNewTimer("ADAPTER::ThermoBaseAlgorithm::SetupThermo");
  Teuchos::TimeMonitor monitor(*t);

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("thermo");

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

//  // get input parameter lists and copy them, because a few parameters are overwritten
  const Teuchos::RCP<Teuchos::ParameterList> ioflags
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->IOParams()));
  const Teuchos::RCP<Teuchos::ParameterList> tdyn
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ThermalDynamicParams()));
//  //const Teuchos::ParameterList& size
//  //  = DRT::Problem::Instance()->ProblemSizeParams();

  // show default parameters of thermo parameter list
  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, *tdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------

  // the default time step size
  tdyn->set<double>("TIMESTEP",prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  tdyn->set<double>("MAXTIME",prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  tdyn->set<int>("NUMSTEP",prbdyn.get<int>("NUMSTEP"));
  // restart
  tdyn->set<int>("RESTARTEVRY",prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  tdyn->set<int>("RESEVRYGLOB",prbdyn.get<int>("UPRES"));
  // get the solver number used for thermal solver
  const int linsolvernumber = tdyn->get<int>("LINEAR_SOLVER");
  // check if the THERMAL solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for thermal solver. Please set LINEAR_SOLVER in THERMAL DYNAMIC to a valid number!");

  // create a linear solver
  Teuchos::RCP<Teuchos::ParameterList> solveparams
    = Teuchos::rcp(new ParameterList());
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // create marching time integrator
  Teuchos::RCP<Thermo> tmpthr;
  switch (timinttype)
  {
  case INPAR::THR::dyna_statics :
  {
    tmpthr = Teuchos::rcp(new THR::TimIntStatics(*ioflags, *tdyn, *xparams,
                                                 actdis, solver, output));
    break;
  }
  case INPAR::THR::dyna_onesteptheta :
  {
    tmpthr = Teuchos::rcp(new THR::TimIntOneStepTheta(*ioflags, *tdyn, *xparams,
                                                      actdis, solver, output));
    break;
  }
  case INPAR::THR::dyna_genalpha :
  {
    dserror("Sorry, GenAlpha not yet available");
    tmpthr = Teuchos::rcp(new THR::TimIntGenAlpha(*ioflags, *tdyn, *xparams,
                                                  actdis, solver, output));
    break;
  }
  case INPAR::THR::dyna_expleuler :
  {
    tmpthr = Teuchos::rcp(new THR::TimIntExplEuler(*ioflags, *tdyn, *xparams,
                                                   actdis, solver, output));
    break;
  }
  default :
    dserror("unknown time integration scheme '%s'", timinttype);
    break;
  }

  // link/store thermal field solver
  thermo_ = tmpthr;

  // see you
  return;
}

/*----------------------------------------------------------------------*
 | Integrate                                                bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::Thermo::Integrate()
{
  // a few parameters
  double time = GetTime();
  const double timeend = GetTimeEnd();
  const double timestepsize = GetTimeStepSize();
  int step = GetTimeStep();
  const int stepend = GetTimeNumStep();

  // loop ahead --- if timestepsize>0
  while ( ((time + (1.e-10)*GetTimeStepSize())< timeend) and (step < stepend) )
  {
    IntegrateStep();

    // update
    Update();
    time +=  timestepsize;
    step += 1;

    // print step summary
    PrintStep();

//    // older version talk to user
//    fprintf(stdout,
//            "Finalised: step %6d"
//            " | nstep %6d"
//            " | time %-14.8E"
//            " | dt %-14.8E\n",
//            step, stepend, time, timestepsize);
//    // print a beautiful line made exactly of 80 dashes
//    fprintf(stdout,
//            "--------------------------------------------------------------"
//            "------------------\n");
//    // do it, print now!
//    fflush(stdout);

    // talk to disk
    Output();
  }

  // print monitoring of time consumption
  Teuchos::TimeMonitor::summarize();

  // Jump you f***ers
  return;
}

/*----------------------------------------------------------------------*/
