/*----------------------------------------------------------------------*/
/*!
\file particle_timint_genalpha.cpp
\brief Implicit particle time integration scheme (backward Euler?)

\level 2

<pre>
\maintainer Alessandro Cattabiani
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_genalpha.H"

#include "particle_algorithm.H"
#include "particleMeshFree_interaction.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "particle_algorithm.H"
#include "particle_heatSource.H"
/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntGenAlpha::TimIntGenAlpha(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntImpl(ioparams, particledynparams, xparams, actdis, output)
{

  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntGenAlpha::Init()
{
  // check
  if (particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::MeshFree)
  {
    dserror("this interaction model is not combined with the hybrid divergence free scheme");
  }

  // call base class init
  PARTICLE::TimIntImpl::Init();

  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, restDensity_));



}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntGenAlpha::IntegrateStep()
{


  return 0;
}

/*----------------------------------------------------------------------*/
// overload to determine the initial accelerations
void PARTICLE::TimIntGenAlpha::DetermineMassDampConsistAccel()
{
  interHandler_->Init(disn_, veln_, radiusn_, mass_, densityn_,specEnthalpyn_,pressure_, temperature_, densityapproxn_,  stepn_);
  interHandler_->Inter_pvp_acc(accn_);

  return;
}

