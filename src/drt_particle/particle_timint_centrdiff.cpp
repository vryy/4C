/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.cpp
\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

\level 2

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_centrdiff.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "particle_utils.H"
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
PARTICLE::TimIntCentrDiff::TimIntCentrDiff(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntExpl
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  )
{
  // DetermineMassDampConsistAccel() is called at the end of Algorithm::Init() after proper setup of the problem
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntCentrDiff::Init()
{

  // call base class init
  PARTICLE::TimIntExpl::Init();

  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();

  switch(particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree:
  {
    interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, restDensity_));
    break;
  }
  case INPAR::PARTICLE::Normal_DEM:
  case INPAR::PARTICLE::Normal_DEM_thermo:
  case INPAR::PARTICLE::NormalAndTang_DEM:
  {
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerDEM(discret_, particle_algorithm_, particleparams));
    break;
  }
  case INPAR::PARTICLE::Normal_MD:
    dserror("central difference time integrator cannot be combined with molecular dynamics collision mechanism");
  break;
  default:
  {
    if(myrank_ == 0)
      std::cout << "central difference time integrator is not combined with a collision model" << std::endl;
    break;
  }
  }

  // check for validity of input data
  if(collhandler_ != Teuchos::null)
  {
    const int numparticle = (*radius_)(0)->GlobalLength();
    if(numparticle > 0)
    {
      double maxradius = 1.0e12;
      (*radius_)(0)->MaxValue(&maxradius);
      if(maxradius > collhandler_->GetMaxRadius())
        dserror("Input parameter MAX_RADIUS invalid");

      double minradius = -1.0;
      (*radius_)(0)->MinValue(&minradius);
      if(minradius < collhandler_->GetMinRadius())
        dserror("Input parameter MIN_RADIUS invalid");
    }
  }

  // simple check if the expansion speed is too elevated
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  if (extParticleMat != NULL)
  {
    // extract the interesting values
    const double min_density = extParticleMat->initDensity_;//for now particle densities at the beginning are all equal. Since dismembering increases density it is also equal to the min_density. It can change in the future tho
    const double CPS = extParticleMat->CPS_;
    const double inv_CPS = 1/CPS;
    const double CPL = extParticleMat->CPL_;
    const double inv_CPL = 1/CPL;
    const double thermalExpansionS = extParticleMat->thermalExpansionS_;
    const double thermalExpansionL = extParticleMat->thermalExpansionL_;
    const double thermalExpansionT = extParticleMat->thermalExpansionT_;

    // extract the boundaries
    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    const double v_max = particleparams.get<double>("MAX_VELOCITY");
    const double r_max = particleparams.get<double>("MAX_RADIUS");

    // loop over the heat sources
    const std::list <Teuchos::RCP<HeatSource> > heatSources = particle_algorithm_->HeatSources();
    double max_QDot = 0;
    for (std::list <Teuchos::RCP<HeatSource> >::const_iterator iHS = heatSources.begin(); iHS != heatSources.end(); ++iHS)
    {
      if ((*iHS)->QDot_> max_QDot)
        max_QDot += std::abs((*iHS)->QDot_);
    }
    const double deltaSpecEnthalpy = max_QDot * (*dt_)[0]/min_density;

    if(r_max<0.0)
      dserror("A positive value of r_max is required!");
    double volume = PARTICLE::Utils::Radius2Volume(r_max);
    volume *= inv_CPS * thermalExpansionS * deltaSpecEnthalpy + 1;
    double v_max_thermo = 2*(PARTICLE::Utils::Volume2Radius(volume) - r_max);

    volume = PARTICLE::Utils::Radius2Volume(r_max);
    volume *= thermalExpansionT * deltaSpecEnthalpy + 1;
    v_max_thermo = std::max(v_max_thermo,2*(PARTICLE::Utils::Volume2Radius(volume) - r_max));

    volume = PARTICLE::Utils::Radius2Volume(r_max);
    volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpy + 1;
    v_max_thermo = std::max(v_max_thermo,2*(PARTICLE::Utils::Volume2Radius(volume) - r_max));

    if (v_max_thermo > 0.01 * v_max)
      dserror("WARNING! The expansion speed of the particle radii is bigger that 1\% of MAX_VELOCITY. Dismembered particles can explode");
  }

  /////////////
  // fast check of the thermodynamic heat exchange
  //    std::cout << "cheap temperature change to test\n";
  //    (*specEnthalpyn_)[0] += 1e3;
  //    (*(*specEnthalpy_)(0))[0] += 1e3;
  //    std::cout << *specEnthalpyn_ << std::endl;
  //    std::cin.get();
  /////////////
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
    const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
    const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

    // new velocities \f$V_{n+1/2}\f$
    //TODO Christoph: check, if the "postponed" subtraction of the accn_ term causes problems for SPH
    // In DEM, this procedure is known to reduce the temporal convergence order if velocity-propoertional damping terms exist.
    veln_->Update(dthalf, *(*acc_)(0), 1.0);

    // new displacements \f$D_{n+1}\f$
    disn_->Update(dt, *veln_, 1.0);

    // apply Dirichlet BCs
    ApplyDirichletBC(timen_, disn_, Teuchos::null, Teuchos::null, false);
    ApplyDirichletBC(timen_-dthalf, Teuchos::null, veln_, Teuchos::null, false);

    // define vector for contact force and moment
    Teuchos::RCP<Epetra_Vector> f_contact = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> m_contact = Teuchos::null;

    // total internal energy (elastic spring and potential energy)
    intergy_ = 0.0;
    //---------------------Compute Collisions-----------------------
    if(collhandler_ != Teuchos::null)
    {
      // new angular-velocities \f$ang_V_{n+1/2}\f$
      angVeln_->Update(dthalf, *(*angAcc_)(0), 1.0);

      // initialize vectors for contact force and moment
      f_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);
      m_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);

      if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo)
      {
        collhandler_->Init(disn_, veln_, angVeln_, radiusn_, mass_, densityn_, specEnthalpyn_);
      }
      else
      {
        collhandler_->Init(disn_, veln_, angVeln_, (*radius_)(0), mass_);
      }

      intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact, specEnthalpyDotn_);
    }
    //--------------------------------------------------------------

    if (interHandler_ != Teuchos::null)
    {
      // the density update scheme is equal to the acceleration update scheme. It can change at your will

      interHandler_->Init(stepn_, disn_, veln_, radiusn_, mass_, specEnthalpyn_, temperature_, densityn_, pressure_);

      // direct update of the accelerations
      // mass balance
      interHandler_->Inter_pvp_densityDot(densityDotn_);
      interHandler_->Inter_pvw_densityDot(densityDotn_);
      // momentum balance
      interHandler_->Inter_pvp_acc(accn_);
      interHandler_->Inter_pvw_acc(accn_);
      // heat balance
      interHandler_->Inter_pvp_specEnthalpyDot(specEnthalpyDotn_);
      interHandler_->Inter_pvhs_specEnthalpyDot(specEnthalpyDotn_);

      // do we need the second round and the surface tension? here we decide
      if (particle_algorithm_->ExtParticleMat()->surfaceVoidTension_ != 0)
      {
        Teuchos::RCP<Epetra_Vector> colorFieldGradientn = LINALG::CreateVector(*discret_->DofRowMap(),true);
        // compute the color field gradient
        interHandler_->Inter_pvp_colorFieldGradient(colorFieldGradientn);
        // update the ParticleMeshFreeData
        interHandler_->SetStateVector(colorFieldGradientn, ColorFieldGradient);
        // compute the surface tension
        interHandler_->Inter_pvp_surfaceTensionCFG(accn_);
      }

      // clear vectors, keep memory
      interHandler_->Clear();
    }
    else
      ComputeAcc(f_contact, m_contact, accn_, angAccn_);

    //--- update with the new accelerations ---//

    // update of end-velocities \f$V_{n+1}\f$
    veln_->Update(dthalf, *accn_, 1.0);

    switch (particle_algorithm_->ParticleInteractionType())
    {
    case INPAR::PARTICLE::MeshFree :
    {
      //Explicit Euler for Density
      densityn_->Update(dt, *densityDotn_, 1.0);
    }// no break
    case INPAR::PARTICLE::Normal_DEM_thermo :
    {
      specEnthalpyn_->Update(dt, *specEnthalpyDotn_, 1.0);
      break;
    }
    default :
      break;
    }


    if(collhandler_ != Teuchos::null)
    {
      angVeln_->Update(dthalf,*angAccn_,1.0);
      // for visualization of orientation vector
      if(writeorientation_)
        RotateOrientVector(dt);
    }

    // apply Dirichlet BCs
    ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

    return 0;
}
