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
#include "particle_contact.H"
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
  INPAR::PARTICLE::ContactStrategy contact_strategy =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");

  switch(contact_strategy)
  {
  case INPAR::PARTICLE::Normal_DEM:
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
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  // Temperatures are slaves of displacements meaning that there is a time-lag
  // of one time step because ComputeTemperatures is based on disn_ {n+1}.
  // It can change in the future, just invertet their order here

  ComputeDisplacements();

  if (trg_temperature_)
    ComputeTemperatures();

  return 0;
}

/*----------------------------------------------------------------------*/
/* integrate step - temperature related physics*/
int PARTICLE::TimIntCentrDiff::ComputeTemperatures()
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
      for(int i_particle=0; i_particle < currbin->NumNode(); i_particle++)
      {
        // position of the current particle
        const DRT::Node* activeNode = particles[i_particle];

        const int gidDof = discret_->Dof(activeNode,0);
        const int lidDof = disn_->Map().LID(gidDof);

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

          // update temperature
          double temp_LH_increase = (HSQDot * dt)/(*mass_)[lidNode];
          if ((*SL_latent_heat_)[lidNode] == 0)
          {
            double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPS;
            if (newTemperatureNoPhaseTransition <= SL_transitionTemperature_)
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            else
            {
              double deltaLatentHeat = temp_LH_increase-(SL_transitionTemperature_-(*temperaturen_)[lidNode])*CPS_;
              if (deltaLatentHeat <= SL_latent_heat_max_)
              {
                (*SL_latent_heat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = SL_transitionTemperature_;
              }
              else
              {
                (*temperaturen_)[lidNode] = SL_transitionTemperature_ + (deltaLatentHeat - SL_latent_heat_max_)*inv_CPL;
                (*SL_latent_heat_)[lidNode] = SL_latent_heat_max_;
              }
            }
          }
          else if ((*SL_latent_heat_)[lidNode] == SL_latent_heat_max_)
          {
            double newTemperatureNoPhaseTransition = (*temperaturen_)[lidNode] + temp_LH_increase*inv_CPL;
            if (newTemperatureNoPhaseTransition >= SL_transitionTemperature_)
              (*temperaturen_)[lidNode] = newTemperatureNoPhaseTransition;
            else
            {
              double deltaLatentHeat = temp_LH_increase-(SL_transitionTemperature_-(*temperaturen_)[lidNode])*CPL_;
              if (-deltaLatentHeat <= SL_latent_heat_max_)
              {
                (*SL_latent_heat_)[lidNode] += deltaLatentHeat;
                (*temperaturen_)[lidNode] = SL_transitionTemperature_;
              }
              else
              {
                (*temperaturen_)[lidNode] = SL_transitionTemperature_ + (deltaLatentHeat + SL_latent_heat_max_)*inv_CPS;
                (*SL_latent_heat_)[lidNode] = 0;
              }
            }
          }
          else if ((*SL_latent_heat_)[lidNode] < SL_latent_heat_max_ && (*SL_latent_heat_)[lidNode] > 0)
          {
            double newLatentHeatKeepTransitioning = (*SL_latent_heat_)[lidNode] + temp_LH_increase;
            if (newLatentHeatKeepTransitioning > SL_latent_heat_max_)
            {
              (*temperaturen_)[lidNode] = SL_transitionTemperature_ + (newLatentHeatKeepTransitioning - SL_latent_heat_max_)*inv_CPL;
              (*SL_latent_heat_)[lidNode] = SL_latent_heat_max_;
            }
            else if (newLatentHeatKeepTransitioning < 0)
            {
              (*temperaturen_)[lidNode] = SL_transitionTemperature_ + newLatentHeatKeepTransitioning*inv_CPS;
              (*SL_latent_heat_)[lidNode] = 0;
            }
            else
              (*SL_latent_heat_)[lidNode] = newLatentHeatKeepTransitioning;
          }
          else
            dserror("latent heat %d out of range (0 : %d). There is a nasty bug in the code?",(*SL_latent_heat_)[lidNode],SL_latent_heat_max_);
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
