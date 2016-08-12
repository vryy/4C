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
#include "../drt_mat/particleAMmat.H"
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
  case INPAR::PARTICLE::Normal_DEM:
  case INPAR::PARTICLE::Normal_DEM_thermo:
  case INPAR::PARTICLE::NormalAndTang_DEM:
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerDEM(discret_, particle_algorithm_, particleparams));
  break;
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
    const int numparticle = radius_->GlobalLength();
    if(numparticle > 0)
    {
      double maxradius = 1.0e12;
      radius_->MaxValue(&maxradius);
      if(maxradius > collhandler_->GetMaxRadius())
        dserror("Input parameter MAX_RADIUS invalid");

      double minradius = -1.0;
      radius_->MinValue(&minradius);
      if(minradius < collhandler_->GetMinRadius())
        dserror("Input parameter MIN_RADIUS invalid");
    }
  }

  // simple check if the expansion speed is too elevated
  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo)
  {
    const int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particleAMmat);
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleAMmat* actmat_derived = static_cast<const MAT::PAR::ParticleAMmat*>(mat);

    const std::map<int,Teuchos::RCP<HeatSource> > heatSources = particle_algorithm_->HeatSources();
    double max_HSQDot = 0;
    for (std::map<int,Teuchos::RCP<HeatSource> >::const_iterator iHS = heatSources.begin(); iHS != heatSources.end(); ++iHS)
    {
      if (iHS->second->HSQDot_> max_HSQDot)
        max_HSQDot = iHS->second->HSQDot_;
    }
    const double min_density = actmat_derived->density_;//for now particle densities at the beginning are all equal. Since dismembering increases density it is also equal to the min_density. It can change in the future tho

    const double delta_SL = max_HSQDot * (*dt_)[0]/min_density;
    const double delta_S = delta_SL/CPS_;
    const double delta_L = delta_SL/CPL_;

    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    const double v_max = particleparams.get<double>("MAX_VELOCITY");
    const double r_max = particleparams.get<double>("MAX_RADIUS");
    double r_max_new = r_max;
    RadiusUpdater(S_thermalExpansion_,delta_S,r_max_new);
    double v_max_thermo = 2*(r_max_new - r_max);
    r_max_new = r_max;
    RadiusUpdater(SL_thermalExpansion_,delta_SL,r_max_new);
    v_max_thermo = std::max(v_max_thermo,2*(r_max_new - r_max));
    r_max_new = r_max;
    RadiusUpdater(SL_thermalExpansion_,delta_L,r_max_new);
    v_max_thermo = std::max(v_max_thermo,2*(r_max_new - r_max));

    if (v_max_thermo > 0.01 * v_max)
      dserror("WARNING! The expansion speed of the particle radii is bigger that 1\% of MAX_VELOCITY. Dismembered particles can explode");
  }
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

  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo)
  {
    ComputeThermodynamics();
  }

  return 0;
}

/*----------------------------------------------------------------------*/
/* integrate step - temperature related physics*/
int PARTICLE::TimIntCentrDiff::ComputeThermodynamics()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$

  // set up the new temperature vector
  temperaturen_->Update(1.0, *(*temperature_)(0), 0.0);

  // extract the map
  const std::map<int, std::list<Teuchos::RCP<HeatSource> > >& bins2heatSources = particle_algorithm_->GetBins2HeatSources();

  double inv_CPS = 1/CPS_;
  double inv_CPL = 1/CPL_;

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
      const double* HSZone_minVer = &(*i_HS)->HSZone_minVer_[0];
      const double* HSZone_maxVer = &(*i_HS)->HSZone_maxVer_[0];
      const double HSQDot = (*i_HS)->HSQDot_;

      // cycle over the nodes of the bin
      const int N_particlesNoAdditions = currbin->NumNode();
      // set the max of the cycle to avoid the added particles
      for(int i_particle=0; i_particle < N_particlesNoAdditions; i_particle++)
      {

        // position of the current particle
        DRT::Node* activeNode = particles[i_particle];

        int gidDof = discret_->Dof(activeNode,0);
        int lidDof = disn_->Map().LID(gidDof);

        if (HSZone_minVer[0]<=(*disn_)[lidDof  ] &&
            HSZone_minVer[1]<=(*disn_)[lidDof+1] &&
            HSZone_minVer[2]<=(*disn_)[lidDof+2] &&
            HSZone_maxVer[0]>=(*disn_)[lidDof  ] &&
            HSZone_maxVer[1]>=(*disn_)[lidDof+1] &&
            HSZone_maxVer[2]>=(*disn_)[lidDof+2])
        {
          // find the node position
          const int gidNode = activeNode->Id();
          const int lidNode = temperaturen_->Map().LID(gidNode);
          if(lidNode < 0)
            dserror("lidNode is not on this proc ");

          // update temperatures, SL_latent_heat, densities, and radii
          const double temp_LH_increase = (HSQDot * dt)/(*density_)[lidNode];

          // solid state start
          if ((*SL_latent_heat_)[lidNode] == 0)
          {
            const double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPS;

            // no phase transition
            if (newTemperatureNoPhaseTransition <= SL_transitionTemperature_)
            {
              // delta temperature NO transition
              const double S_deltaTemperature = newTemperatureNoPhaseTransition - (*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(S_thermalExpansion_, S_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(S_thermalExpansion_, S_deltaTemperature, (*radius_)[lidNode]);
              // temperature update
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            }

            // a (wild) phase transition appears! (solid -> liquid)
            else
            {
              // delta temperature BEFORE phase transition
              const double S_deltaTemperature = SL_transitionTemperature_-(*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(S_thermalExpansion_, S_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(S_thermalExpansion_, S_deltaTemperature, (*radius_)[lidNode]);
              // delta latent heat DURING phase transition
              const double deltaLatentHeat = temp_LH_increase - S_deltaTemperature * CPS_;

              // still in the phase transition after this the time step
              if (deltaLatentHeat <= SL_latent_heat_max_)
              {
                // density and radius update
                DensityUpdater(SL_thermalExpansion_, deltaLatentHeat, (*density_)[lidNode]);
                RadiusUpdater(SL_thermalExpansion_, deltaLatentHeat, (*radius_)[lidNode]);
                // temperature and latent heat updates
                (*SL_latent_heat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = SL_transitionTemperature_;
              }

              // rare case where the phase transition happens in 1 time step
              else
              {
                // density and radius update
                DensityUpdater(SL_thermalExpansion_, SL_latent_heat_max_, (*density_)[lidNode]);
                RadiusUpdater(SL_thermalExpansion_, SL_latent_heat_max_, (*radius_)[lidNode]);
                // delta temperature AFTER phase transition
                const double L_deltaTemperature = (deltaLatentHeat - SL_latent_heat_max_)*inv_CPL;
                // density and radius update
                DensityUpdater(L_thermalExpansion_, L_deltaTemperature, (*density_)[lidNode]);
                RadiusUpdater(L_thermalExpansion_, L_deltaTemperature, (*radius_)[lidNode]);
                // temperature and latent heat updates
                (*temperaturen_)[lidNode] = SL_transitionTemperature_ + L_deltaTemperature;
                (*SL_latent_heat_)[lidNode] = SL_latent_heat_max_;
              }
            }
          }

          //liquid state start
          else if ((*SL_latent_heat_)[lidNode] == SL_latent_heat_max_)
          {
            const double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPL;

            // no phase transition
            if (newTemperatureNoPhaseTransition >= SL_transitionTemperature_)
            {
              // delta temperature NO transition
              const double L_deltaTemperature = newTemperatureNoPhaseTransition - (*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(L_thermalExpansion_, L_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(L_thermalExpansion_, L_deltaTemperature, (*radius_)[lidNode]);
              // temperature update
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            }

            // a (wild) phase transition appears! (liquid -> solid)
            else
            {
              // delta temperature BEFORE phase transition
              const double L_deltaTemperature = SL_transitionTemperature_-(*temperaturen_)[lidNode];
              // density and radius update
              DensityUpdater(L_thermalExpansion_, L_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(L_thermalExpansion_, L_deltaTemperature, (*radius_)[lidNode]);
              // delta latent heat DURING phase transition
              const double deltaLatentHeat = temp_LH_increase-(SL_transitionTemperature_-(*temperaturen_)[lidNode])*CPL_;

              // still in the phase transition after this the time step
              if (-deltaLatentHeat <= SL_latent_heat_max_)
              {
                // density and radius update
                DensityUpdater(SL_thermalExpansion_, deltaLatentHeat, (*density_)[lidNode]);
                RadiusUpdater(SL_thermalExpansion_, deltaLatentHeat, (*radius_)[lidNode]);
                // temperature and latent heat updates
                (*SL_latent_heat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = SL_transitionTemperature_;
              }

              // rare case where the phase transition happens in 1 time step
              else
              {
                // density and radius update
                DensityUpdater(SL_thermalExpansion_, SL_latent_heat_max_, (*density_)[lidNode]);
                RadiusUpdater(SL_thermalExpansion_, SL_latent_heat_max_, (*radius_)[lidNode]);
                // delta temperature AFTER phase transition
                const double S_deltaTemperature = (deltaLatentHeat + SL_latent_heat_max_)*inv_CPS;
                // density and radius update
                DensityUpdater(L_thermalExpansion_, L_deltaTemperature, (*density_)[lidNode]);
                RadiusUpdater(L_thermalExpansion_, L_deltaTemperature, (*radius_)[lidNode]);
                // temperature and latent heat updates
                (*temperaturen_)[lidNode] = SL_transitionTemperature_ + S_deltaTemperature;
                (*SL_latent_heat_)[lidNode] = 0;
              }
            }
          }

          // transition start (solid <-> liquid)
          else if ((*SL_latent_heat_)[lidNode] < SL_latent_heat_max_ && (*SL_latent_heat_)[lidNode] > 0)
          {
            const double newLatentHeatKeepTransitioning = (*SL_latent_heat_)[lidNode] + temp_LH_increase;

            // transition -> liquid
            if (newLatentHeatKeepTransitioning > SL_latent_heat_max_)
            {
              // delta temperature DURING phase transition
              const double deltaLatentHeat = SL_latent_heat_max_ - (*SL_latent_heat_)[lidNode];
              // density and radius update
              DensityUpdater(SL_thermalExpansion_, deltaLatentHeat, (*density_)[lidNode]);
              RadiusUpdater(SL_thermalExpansion_, deltaLatentHeat, (*radius_)[lidNode]);
              // delta latent heat AFTER phase transition
              const double L_deltaTemperature = (newLatentHeatKeepTransitioning - SL_latent_heat_max_) * inv_CPL;
              // density and radius update
              DensityUpdater(L_thermalExpansion_, L_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(L_thermalExpansion_, L_deltaTemperature, (*radius_)[lidNode]);
              // temperature and latent heat updates
              (*temperaturen_)[lidNode] = SL_transitionTemperature_ + L_deltaTemperature;
              (*SL_latent_heat_)[lidNode] = SL_latent_heat_max_;
            }

            // transition -> solid
            else if (newLatentHeatKeepTransitioning < 0)
            {
              // delta temperature DURING phase transition
              const double deltaLatentHeat = -(*SL_latent_heat_)[lidNode];
              // density and radius update
              DensityUpdater(SL_thermalExpansion_, deltaLatentHeat, (*density_)[lidNode]);
              RadiusUpdater(SL_thermalExpansion_, deltaLatentHeat, (*radius_)[lidNode]);
              // delta latent heat AFTER phase transition
              const double S_deltaTemperature = newLatentHeatKeepTransitioning * inv_CPS;
              // density and radius update
              DensityUpdater(S_thermalExpansion_, S_deltaTemperature, (*density_)[lidNode]);
              RadiusUpdater(S_thermalExpansion_, S_deltaTemperature, (*radius_)[lidNode]);
              // temperature and latent heat updates
              (*temperaturen_)[lidNode] = SL_transitionTemperature_ + S_deltaTemperature;
              (*SL_latent_heat_)[lidNode] = 0;
            }

            // we stay in the transition zone
            else
            {
              // delta temperature DURING phase transition (this is kept just for the sake of readability)
              const double deltaLatentHeat = temp_LH_increase;
              // density and radius update
              DensityUpdater(SL_thermalExpansion_, deltaLatentHeat, (*density_)[lidNode]);
              RadiusUpdater(SL_thermalExpansion_, deltaLatentHeat, (*radius_)[lidNode]);
              // latent heat update
              (*SL_latent_heat_)[lidNode] = newLatentHeatKeepTransitioning;
            }
          }

          // other cases, did I forget anything?
          else
          {
            dserror("latent heat %d out of range (0 : %d). There is a nasty bug in the code?",(*SL_latent_heat_)[lidNode],SL_latent_heat_max_);
          }
        }
      }
    }
  }
  return 0;
}

/*----------------------------------------------------------------------*/
/* integrate step - displacement related physics*/
int PARTICLE::TimIntCentrDiff::ComputeDisplacements()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
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
  ang_veln_->Update(1.0, *(*ang_vel_)(0), 0.0);
  ang_veln_->Update(dthalf, *(*ang_acc_)(0), 1.0);

  // initialize vectors for contact force and moment
  f_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  m_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);

  SetStatesForCollision();

  intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact);
  }
  //--------------------------------------------------------------

  ComputeAcc(f_contact, m_contact, accn_, ang_accn_);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);
  if(collhandler_ != Teuchos::null)
  {
    ang_veln_->Update(dthalf,*ang_accn_,1.0);
    // for visualization of orientation vector
    if(writeorientation_)
      RotateOrientVector(dt);
  }

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  return 0;
}


