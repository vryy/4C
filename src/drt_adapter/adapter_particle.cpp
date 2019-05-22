/*----------------------------------------------------------------------*/
/*!

\brief Particle field adapter

\level 1

\maintainer  Sebastian Fuchs


*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/13 |
 *----------------------------------------------------------------------*/
#include "adapter_particle.H"
#include "../drt_particle/particle_timint_centrdiff.H"
#include "../drt_particle/particle_timint_expleuler.H"
#include "../drt_particle/particle_timint_rk.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 |                                                          ghamm 12/13 |
 *----------------------------------------------------------------------*/
ADAPTER::ParticleBaseAlgorithm::ParticleBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  SetupTimInt(prbdyn, actdis);
}


/*----------------------------------------------------------------------*
 | setup of particle time integration                       ghamm 12/13 |
 *----------------------------------------------------------------------*/
void ADAPTER::ParticleBaseAlgorithm::SetupTimInt(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ADAPTER::ParticleBaseAlgorithm::SetupParticle");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();

  // get input parameter lists and copy them, because a few parameters are overwritten
  const Teuchos::ParameterList& ioflags = DRT::Problem::Instance()->IOParams();
  const Teuchos::RCP<Teuchos::ParameterList> partdyn =
      Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ParticleParamsOld()));

  // show default parameters of particle parameter list
  if ((actdis->Comm()).MyPID() == 0) DRT::INPUT::PrintDefaultParameters(IO::cout, *partdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  partdyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  partdyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  partdyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  partdyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  // switch to different time integrators
  INPAR::PARTICLEOLD::DynamicType timinttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::DynamicType>(*partdyn, "DYNAMICTYP");

  // create marching time integrator
  Teuchos::RCP<Particle> tmppart;
  switch (timinttype)
  {
    case INPAR::PARTICLEOLD::dyna_expleuler:
    {
      tmppart =
          Teuchos::rcp(new PARTICLE::TimIntExplEuler(ioflags, *partdyn, *xparams, actdis, output));
      break;
    }
    case INPAR::PARTICLEOLD::dyna_centrdiff:
    {
      tmppart =
          Teuchos::rcp(new PARTICLE::TimIntCentrDiff(ioflags, *partdyn, *xparams, actdis, output));
      break;
    }
    case INPAR::PARTICLEOLD::dyna_rk2:
    case INPAR::PARTICLEOLD::dyna_rk4:
    {
      tmppart = Teuchos::rcp(new PARTICLE::TimIntRK(ioflags, *partdyn, *xparams, actdis, output));
      break;
    }
    default:
    {
      dserror("unknown time integration scheme '%s'", timinttype);
      break;
    }
  }

  // store particle field
  particle_ = tmppart;

  return;
}  // SetupTimInt()
