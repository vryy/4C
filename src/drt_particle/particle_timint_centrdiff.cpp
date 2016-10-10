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
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "particle_algorithm.H"
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
  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();

  switch(particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree:
  {
    //interhandler_ = Teuchos::rcp(new PARTICLE::MeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams));
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
  }
  break;
  }

  // call base class init
  PARTICLE::TimIntExpl::Init();

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
  std::cout << "Warning! You have erased the expansion check!\n";
  /*
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  if (extParticleMat != NULL)
  {
    const std::map<int,Teuchos::RCP<HeatSource> > heatSources = particle_algorithm_->HeatSources();
    double max_QDot = 0;
    for (std::map<int,Teuchos::RCP<HeatSource> >::const_iterator iHS = heatSources.begin(); iHS != heatSources.end(); ++iHS)
    {
      if (iHS->second->QDot_> max_QDot)
        max_QDot = iHS->second->QDot_;
    }
    const double min_density = extParticleMat->initDensity_;//for now particle densities at the beginning are all equal. Since dismembering increases density it is also equal to the min_density. It can change in the future tho

    const double delta_SL = max_QDot * (*dt_)[0]/min_density;
    const double delta_S = delta_SL/extParticleMat->CPS_;
    const double delta_L = delta_SL/extParticleMat->CPL_;

    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    const double v_max = particleparams.get<double>("MAX_VELOCITY");
    const double r_max = particleparams.get<double>("MAX_RADIUS");
    double r_max_new = r_max;
    RadiusUpdater(extParticleMat->thermalExpansionS_,delta_S,r_max_new);
    double v_max_thermo = 2*(r_max_new - r_max);
    r_max_new = r_max;
    RadiusUpdater(extParticleMat->thermalExpansionSL_,delta_SL,r_max_new);
    v_max_thermo = std::max(v_max_thermo,2*(r_max_new - r_max));
    r_max_new = r_max;
    RadiusUpdater(extParticleMat->thermalExpansionSL_,delta_L,r_max_new);
    v_max_thermo = std::max(v_max_thermo,2*(r_max_new - r_max));

    if (v_max_thermo > 0.01 * v_max)
      dserror("WARNING! The expansion speed of the particle radii is bigger that 1\% of MAX_VELOCITY. Dismembered particles can explode");
  }
  */
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  // Temperatures are slaves of displacements meaning that there is a time-lag
  // of one time step because ComputeTemperatures is based on disn_ {n+1}.
  // It can change in the future, just invert their order here
  ComputeDisplacements();

  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo ||
      particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::MeshFree)
    ComputeThermodynamics();

  return 0;
}

/*----------------------------------------------------------------------*/
/* integrate step - thermodynamics*/
void PARTICLE::TimIntCentrDiff::ComputeThermodynamics1()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$

  // extract the material properties
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  const double specEnthalpyST = extParticleMat->SpecEnthalpyST();
  const double specEnthalpyTL = extParticleMat->SpecEnthalpyTL();
  const double CPS = extParticleMat->CPS_;
  const double inv_CPS = 1/CPS;
  const double CPL = extParticleMat->CPL_;
  const double inv_CPL = 1/CPL;
  const double latentHeatSL = extParticleMat->latentHeatSL_;
  const double thermalExpansionS = extParticleMat->thermalExpansionS_;
  const double thermalExpansionL = extParticleMat->thermalExpansionL_;
  const double thermalExpansionSL = extParticleMat->thermalExpansionSL_;

  // extract the map
  const std::map<int, std::list<Teuchos::RCP<HeatSource> > >& bins2heatSources = particle_algorithm_->Bins2HeatSources();

  // cycle over the map
  std::map<int, std::list<Teuchos::RCP<HeatSource> > >::const_iterator i_bin;
  for (i_bin = bins2heatSources.begin(); i_bin != bins2heatSources.end(); ++i_bin)
  {
    // extract the current bin
    DRT::Element* currbin = discret_->gElement(i_bin->first);

    if(currbin->Owner() != myrank_)
      dserror("trying to eval a ghost bin :-(");

    // find the nodes of the bin
    DRT::Node** particles = currbin->Nodes();

    // cycle over the heat sources
    std::list<Teuchos::RCP<HeatSource> >::const_iterator i_HS;
    for (i_HS = i_bin->second.begin(); i_HS != i_bin->second.end(); i_HS++)
    {
      // extract vertexes
      const double* minVerZone = &(*i_HS)->minVerZone_[0];
      const double* maxVerZone = &(*i_HS)->maxVerZone_[0];
      const double QDot = (*i_HS)->QDot_;

      // cycle over the nodes of the bin
      const int numNode = currbin->NumNode();
      for(int i_particle=0; i_particle < numNode; i_particle++)
      {
        // position of the current particle
        DRT::Node* activeNode = particles[i_particle];

        int gidDof = discret_->Dof(activeNode,0);
        int lidDof = disn_->Map().LID(gidDof);

        if (minVerZone[0]<=(*disn_)[lidDof  ] &&
            minVerZone[1]<=(*disn_)[lidDof+1] &&
            minVerZone[2]<=(*disn_)[lidDof+2] &&
            maxVerZone[0]>=(*disn_)[lidDof  ] &&
            maxVerZone[1]>=(*disn_)[lidDof+1] &&
            maxVerZone[2]>=(*disn_)[lidDof+2])
        {
          // find the node position
          const int gidNode = activeNode->Id();
          const int lidNode = discret_->NodeRowMap()->LID(gidNode);
          if(lidNode < 0)
            dserror("lidNode is not on this proc ");
          // update specEnthalpy
          (*specEnthalpyn_)[lidNode] = (QDot * dt)/(*densityn_)[lidNode];
        }
      }
    }
  }

  // update the other state vectors (\rho and R)
  for (int lidNode = 0; lidNode < discret_->NumMyRowNodes(); ++lidNode)
  {
    const double oldSpecEnthalpy = (*(*specEnthalpy_)(0))[lidNode];
    const double newSpecEnthalpy = (*specEnthalpyn_)[lidNode];
    // skip in case the specEnthalpy did not change
    if (newSpecEnthalpy != oldSpecEnthalpy)
    {
      // compute the current volume
      double volume = Radius2Volume((*radiusn_)[lidNode]);
      // specEnthalpy difference
      double deltaSpecEnthalpy = newSpecEnthalpy - oldSpecEnthalpy;

      // --- compute the new volume --- //

      // WAS it solid?
      if (oldSpecEnthalpy <= specEnthalpyST)
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
          volume *= EffExpCoeff(inv_CPS * thermalExpansionS,deltaSpecEnthalpy);
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in solid state
          volume *= EffExpCoeff(inv_CPS * thermalExpansionS,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= EffExpCoeff(thermalExpansionSL,latentHeatSL);
          deltaSpecEnthalpy -= latentHeatSL;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpy);
        }
        // it IS transition state
        else
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in solid state
          volume *= EffExpCoeff(inv_CPS * thermalExpansionS,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= EffExpCoeff(thermalExpansionSL,deltaSpecEnthalpy);
        }
      }
      // WAS it liquid?
      else if (oldSpecEnthalpy >= specEnthalpyTL)
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in the transition state
          volume *= EffExpCoeff(thermalExpansionSL,-latentHeatSL);
          deltaSpecEnthalpy -= -latentHeatSL;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpy);
        }
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpy);
        // it IS transition state
        else
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= EffExpCoeff(thermalExpansionSL,deltaSpecEnthalpy);
        }
      }
      // it WAS in transition state
      else
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in transition state
          volume *= EffExpCoeff(thermalExpansionSL,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPS * thermalExpansionS,deltaSpecEnthalpy);
        }
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in transition state
          volume *= EffExpCoeff(thermalExpansionSL,deltaSpecEnthalpyUpToTransition);
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in liquid state
          volume *= EffExpCoeff(inv_CPL * thermalExpansionL,deltaSpecEnthalpy);
        }
        // it IS transition state
        else
          volume *= EffExpCoeff(thermalExpansionSL,deltaSpecEnthalpy);
      }

      // --- compute the new volume --- //

      // updates
      (*radiusn_)[lidNode] = Volume2Radius(volume);
      (*densityn_)[lidNode] = (*mass_)[lidNode]/volume;
    }
  }
}

/*----------------------------------------------------------------------*/
/* integrate step - temperature related physics*/
int PARTICLE::TimIntCentrDiff::ComputeThermodynamics()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$

  // extract the material pointer
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  // extract values for simpler access
  const double CPS = extParticleMat->CPS_;
  const double CPL = extParticleMat->CPL_;
  const double transitionTemperatureSL = extParticleMat->transitionTemperatureSL_;
  const double thermalExpansionS = extParticleMat->thermalExpansionS_;
  const double thermalExpansionL = extParticleMat->thermalExpansionL_;
  const double thermalExpansionSL = extParticleMat->thermalExpansionSL_;
  const double latentHeatSL = extParticleMat->latentHeatSL_;

  // extract the map
  const std::map<int, std::list<Teuchos::RCP<HeatSource> > >& bins2heatSources = particle_algorithm_->Bins2HeatSources();



  double inv_CPS = 1/CPS;
  double inv_CPL = 1/CPL;

  // cycle over the map
  std::map<int, std::list<Teuchos::RCP<HeatSource> > >::const_iterator i_bin;
  for (i_bin = bins2heatSources.begin(); i_bin != bins2heatSources.end(); ++i_bin)
  {
    // extract the current bin
    DRT::Element* currbin = discret_->gElement(i_bin->first);

    if(currbin->Owner() != myrank_)
      dserror("trying to eval a ghost bin :-(");

    // find the nodes of the bin
    DRT::Node** particles = currbin->Nodes();

    // cycle over the heat sources
    std::list<Teuchos::RCP<HeatSource> >::const_iterator i_HS;
    for (i_HS = i_bin->second.begin(); i_HS != i_bin->second.end(); i_HS++)
    {
      // extract vertexes
      const double* minVerZone = &(*i_HS)->minVerZone_[0];
      const double* maxVerZone = &(*i_HS)->maxVerZone_[0];
      const double QDot = (*i_HS)->QDot_;

      // cycle over the nodes of the bin
      const int N_particlesNoAdditions = currbin->NumNode();
      // set the max of the cycle to avoid the added particles
      for(int i_particle=0; i_particle < N_particlesNoAdditions; i_particle++)
      {

        // position of the current particle
        DRT::Node* activeNode = particles[i_particle];

        int gidDof = discret_->Dof(activeNode,0);
        int lidDof = disn_->Map().LID(gidDof);

        if (minVerZone[0]<=(*disn_)[lidDof  ] &&
            minVerZone[1]<=(*disn_)[lidDof+1] &&
            minVerZone[2]<=(*disn_)[lidDof+2] &&
            maxVerZone[0]>=(*disn_)[lidDof  ] &&
            maxVerZone[1]>=(*disn_)[lidDof+1] &&
            maxVerZone[2]>=(*disn_)[lidDof+2])
        {
          // find the node position
          const int gidNode = activeNode->Id();
          const int lidNode = temperaturen_->Map().LID(gidNode);
          if(lidNode < 0)
            dserror("lidNode is not on this proc ");

          // update temperatures, SL_latent_heat, densities, and radii
          const double temp_LH_increase = (QDot * dt)/(*densityn_)[lidNode];

          // solid state start
          if ((*latentHeat_)[lidNode] == 0)
          {
            const double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPS;

            // no phase transition
            if (newTemperatureNoPhaseTransition <= transitionTemperatureSL)
            {
              // delta temperature NO transition
              const double S_deltaTemperature = newTemperatureNoPhaseTransition - (*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionS, S_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionS, S_deltaTemperature, (*radiusn_)[lidNode]);
              // temperature update
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            }

            // a (wild) phase transition appears! (solid -> liquid)
            else
            {
              // delta temperature BEFORE phase transition
              const double S_deltaTemperature = transitionTemperatureSL-(*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionS, S_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionS, S_deltaTemperature, (*radiusn_)[lidNode]);
              // delta latent heat DURING phase transition
              const double deltaLatentHeat = temp_LH_increase - S_deltaTemperature * CPS;

              // still in the phase transition after this the time step
              if (deltaLatentHeat <= latentHeatSL)
              {
                // density and radius update
                DensityUpdater(thermalExpansionSL, deltaLatentHeat, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionSL, deltaLatentHeat, (*radiusn_)[lidNode]);
                // temperature and latent heat updates
                (*latentHeat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = transitionTemperatureSL;
              }

              // rare case where the phase transition happens in 1 time step
              else
              {
                // density and radius update
                DensityUpdater(thermalExpansionSL, latentHeatSL, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionSL, latentHeatSL, (*radiusn_)[lidNode]);
                // delta temperature AFTER phase transition
                const double L_deltaTemperature = (deltaLatentHeat - latentHeatSL)*inv_CPL;
                // density and radius update
                DensityUpdater(thermalExpansionL, L_deltaTemperature, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionL, L_deltaTemperature, (*radiusn_)[lidNode]);
                // temperature and latent heat updates
                (*temperaturen_)[lidNode] = transitionTemperatureSL + L_deltaTemperature;
                (*latentHeat_)[lidNode] = latentHeatSL;
              }
            }
          }

          //liquid state start
          else if ((*latentHeat_)[lidNode] == latentHeatSL)
          {
            const double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPL;

            // no phase transition
            if (newTemperatureNoPhaseTransition >= transitionTemperatureSL)
            {
              // delta temperature NO transition
              const double L_deltaTemperature = newTemperatureNoPhaseTransition - (*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionL, L_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionL, L_deltaTemperature, (*radiusn_)[lidNode]);
              // temperature update
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            }

            // a (wild) phase transition appears! (liquid -> solid)
            else
            {
              // delta temperature BEFORE phase transition
              const double L_deltaTemperature = transitionTemperatureSL-(*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionL, L_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionL, L_deltaTemperature, (*radiusn_)[lidNode]);
              // delta latent heat DURING phase transition
              const double deltaLatentHeat = temp_LH_increase-(transitionTemperatureSL-(*temperaturen_)[lidNode])*CPL;

              // still in the phase transition after this the time step
              if (-deltaLatentHeat <= latentHeatSL)
              {
                // density and radius update
                DensityUpdater(thermalExpansionSL, deltaLatentHeat, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionSL, deltaLatentHeat, (*radiusn_)[lidNode]);
                // temperature and latent heat updates
                (*latentHeat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = transitionTemperatureSL;
              }

              // rare case where the phase transition happens in 1 time step
              else
              {
                // density and radius update
                DensityUpdater(thermalExpansionSL, latentHeatSL, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionSL, latentHeatSL, (*radiusn_)[lidNode]);
                // delta temperature AFTER phase transition
                const double S_deltaTemperature = (deltaLatentHeat + latentHeatSL)*inv_CPS;
                // density and radius update
                DensityUpdater(thermalExpansionL, L_deltaTemperature, (*densityn_)[lidNode]);
                RadiusUpdater(thermalExpansionL, L_deltaTemperature, (*radiusn_)[lidNode]);
                // temperature and latent heat updates
                (*temperaturen_)[lidNode] = transitionTemperatureSL + S_deltaTemperature;
                (*latentHeat_)[lidNode] = 0;
              }
            }
          }

          // transition start (solid <-> liquid)
          else if ((*latentHeat_)[lidNode] < latentHeatSL && (*latentHeat_)[lidNode] > 0)
          {
            const double newLatentHeatKeepTransitioning = (*latentHeat_)[lidNode] + temp_LH_increase;

            // transition -> liquid
            if (newLatentHeatKeepTransitioning > latentHeatSL)
            {
              // delta temperature DURING phase transition
              const double deltaLatentHeat = latentHeatSL - (*latentHeat_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionSL, deltaLatentHeat, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionSL, deltaLatentHeat, (*radiusn_)[lidNode]);
              // delta latent heat AFTER phase transition
              const double L_deltaTemperature = (newLatentHeatKeepTransitioning - latentHeatSL) * inv_CPL;
              // density and radius update
              DensityUpdater(thermalExpansionL, L_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionL, L_deltaTemperature, (*radiusn_)[lidNode]);
              // temperature and latent heat updates
              (*temperaturen_)[lidNode] = transitionTemperatureSL + L_deltaTemperature;
              (*latentHeat_)[lidNode] = latentHeatSL;
            }

            // transition -> solid
            else if (newLatentHeatKeepTransitioning < 0)
            {
              // delta temperature DURING phase transition
              const double deltaLatentHeat = -(*latentHeat_)[lidNode];
              // density and radius update
              DensityUpdater(thermalExpansionSL, deltaLatentHeat, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionSL, deltaLatentHeat, (*radiusn_)[lidNode]);
              // delta latent heat AFTER phase transition
              const double S_deltaTemperature = newLatentHeatKeepTransitioning * inv_CPS;
              // density and radius update
              DensityUpdater(thermalExpansionS, S_deltaTemperature, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionS, S_deltaTemperature, (*radiusn_)[lidNode]);
              // temperature and latent heat updates
              (*temperaturen_)[lidNode] = transitionTemperatureSL + S_deltaTemperature;
              (*latentHeat_)[lidNode] = 0;
            }

            // we stay in the transition zone
            else
            {
              // delta temperature DURING phase transition (this is kept just for the sake of readability)
              const double deltaLatentHeat = temp_LH_increase;
              // density and radius update
              DensityUpdater(thermalExpansionSL, deltaLatentHeat, (*densityn_)[lidNode]);
              RadiusUpdater(thermalExpansionSL, deltaLatentHeat, (*radiusn_)[lidNode]);
              // latent heat update
              (*latentHeat_)[lidNode] = newLatentHeatKeepTransitioning;
            }
          }

          // other cases, did I forget anything?
          else
          {
            dserror("latent heat %d out of range (0 : %d). There is a nasty bug in the code?",(*latentHeat_)[lidNode],latentHeatSL);
          }
        }
      }
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*/
/* integrate step - displacement related physics*/
void PARTICLE::TimIntCentrDiff::ComputeDisplacements()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
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

    SetStatesForCollision();

    intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact);
  }
  //--------------------------------------------------------------

  ComputeAcc(f_contact, m_contact, accn_, angAccn_);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);
  if(collhandler_ != Teuchos::null)
  {
    angVeln_->Update(dthalf,*angAccn_,1.0);
    // for visualization of orientation vector
    if(writeorientation_)
      RotateOrientVector(dt);
  }

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);
}



