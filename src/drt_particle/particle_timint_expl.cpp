/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expl.cpp

\brief Particle time integration with explicit time integration

\level 2

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"

#include "particle_contact.H"
#include "particleMeshFree_interaction.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntExpl::TimIntExpl(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimInt
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  )
{
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntExpl::Init()
{
  // call base class init
  TimInt::Init();

  // initialize collision handler for particles
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  switch(particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::Normal_DEM:
  case INPAR::PARTICLE::Normal_DEM_thermo:
  case INPAR::PARTICLE::NormalAndTang_DEM:
  {
    if(DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) < 0)
      collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerDEM(discret_, particle_algorithm_, particleparams));
    else
      collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerDEMEllipsoids(discret_, particle_algorithm_, particleparams));
    break;
  }
  case INPAR::PARTICLE::MeshFree:
  {
    interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams));
    break;
  }
  case INPAR::PARTICLE::Normal_MD:
  {
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerMD(discret_, particle_algorithm_, particleparams));
    break;
  }
  default:
  {
    if(myrank_ == 0)
      std::cout << "explicit time integrator is not combined with a collision model" << std::endl;
    break;
  }
  }

  // check for validity of input data
  if(collhandler_ != Teuchos::null)
  {
    const int numparticle = (*radius_)(0)->GlobalLength();
    if(numparticle > 0)
    {
      double maxradius = 1.e12;
      (*radius_)(0)->MaxValue(&maxradius);
      if(maxradius > collhandler_->GetMaxRadius())
        dserror("Input parameter MAX_RADIUS invalid");

      double minradius = -1.;
      (*radius_)(0)->MinValue(&minradius);
      if(minradius < collhandler_->GetMinRadius())
        dserror("Input parameter MIN_RADIUS invalid");
    }
  }

  return;
}
