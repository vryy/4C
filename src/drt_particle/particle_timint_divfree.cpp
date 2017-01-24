/*----------------------------------------------------------------------*/
/*!
\file particle_timint_divfree.cpp
\brief Hybrid particle time integration with implicit density and divergence corrections (hybrid)

\level 2

<pre>
\maintainer Alessandro Cattabiani
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_divfree.H"
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
PARTICLE::TimIntDivFree::TimIntDivFree(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntHybrid(ioparams, particledynparams, xparams, actdis, output),
  timeStepType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::TimeStepType>(particledynparams,"TIMESTEPTYPE")),
  correctDivergenceToll_(particledynparams.get<double>("CORRECT_DIVERGENCE_TOLL")),
  correctDivergenceIter_(particledynparams.get<int>("CORRECT_DIVERGENCE_ITER")),
  correctDensityToll_(particledynparams.get<double>("CORRECT_DENSITY_TOLL")),
  correctDensityIter_(particledynparams.get<int>("CORRECT_DENSITY_ITER"))
{
  trg_warmStart_ = true;
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntDivFree::Init()
{
  // check
  if (particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::MeshFree)
  {
    dserror("this interaction model is not combined with the hybrid divergence free scheme");
  }

  // call base class init
  PARTICLE::TimIntHybrid::Init();

  // set up the proper interaction handler
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, restDensity_));

}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntDivFree::IntegrateStep()
{
  // call init for a warm start
  if (trg_warmStart_)
  {
    interHandler_->Clear();
    interHandler_->Init(disn_,veln_,radiusn_,mass_, densityn_,specEnthalpyn_,pressure_, temperature_,  stepn_);
    trg_warmStart_ = false;
  }

  // set the time step size
  const double dt = SetDt();

  // set up the interaction handler updating neighbours and state vectors
  // for the first iteration loop it has already been called

  // pvp interactions (no pvp pressures)
  interHandler_->Inter_pvp_acc(accn_, false);

  // pvw interactions
  interHandler_->Inter_pvw_densityDot(densityDotn_);
  densityn_->Update(dt, *densityDotn_, 1.0);
  interHandler_->SetStateVector(densityn_, Density);
  interHandler_->Inter_pvw_acc(accn_);

  // first guess on the velocity field
  veln_->Update(dt, *accn_, 1.0);

  // fulfill density - restDensity = 0;
  CorrectDensityError(dt);

  // update displacements
  disn_->Update(dt, *veln_, 1.0);

  // better erase everithing, here we update
  interHandler_->Clear();

  // distribute to ParticleMeshFreeData and find neighbours
  interHandler_->Init(disn_,veln_,radiusn_,mass_, densityn_,specEnthalpyn_,pressure_, temperature_,  stepn_);

  // correct divergence error
  CorrectDivergenceError(dt);

  // update velocities in ParticleMeshFreeData
  interHandler_->SetStateVector(veln_, Vel);

  // heat interactions
  // pvp
  interHandler_->Inter_pvp_specEnthalpyDot(specEnthalpyDotn_);
  // pvhs
  interHandler_->Inter_pvhs_specEnthalpyDot(specEnthalpyDotn_);
  specEnthalpyn_->Update(dt, *specEnthalpyDotn_, 1.0);

  return 0;
}


/*--------------------------------------------------------------------------*
 | correct the divergence error (divFree scheme)               katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::TimIntDivFree::CorrectDivergenceError(const double dt)
{
  double densityDotErr = 0;
  int ii;
  for (ii = 0; ii < correctDivergenceIter_; ++ii)
  {
    // compute the densityDot (that should converge towards zero)
    interHandler_->Inter_pvp_densityDot(densityDotn_);

    // did we converge? Let's check the average
    densityDotn_->Norm2(&densityDotErr);
    densityDotErr /= densityDotn_->GlobalLength();
    if (densityDotErr <= correctDivergenceToll_)
    {
      break;
    }


    //std::cout << "The density dot error is = " << densityDotErr << std::endl;
    //std::cout << "Iteration " << ii << "/" << correctDivergenceIter_ << std::endl;
    //std::cin.get();

    // update the local ParticleMeshFreeData of the processor
    interHandler_->SetStateVector(densityDotn_, DensityDot);

    // create and compute the divergence free pressure field
    Teuchos::RCP<Epetra_Vector> divFreePressureAcc = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
    interHandler_->Inter_pvp_divFreePressureAcc(divFreePressureAcc, dt);

    // update the velocity
    veln_->Update(dt, *divFreePressureAcc, 1.0);
  }

  std::cout << "Divergence corrected with " << ii << "/" << correctDivergenceIter_ << " iterations\n";
  std::cout << "Final density dot error = " << densityDotErr << std::endl;
}


/*--------------------------------------------------------------------------*
 | correct the density error (divFree scheme)                  katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::TimIntDivFree::CorrectDensityError(const double dt)
{
  double densityErr = 0;
  int ii;
  for (ii = 0; ii < correctDensityIter_; ++ii)
  {
    // compute the densityDot (that should converge towards zero)
    densityn_->PutScalar(0.0);
    interHandler_->MF_mW(densityn_);
    interHandler_->SetStateVector(densityn_, Density);

    Teuchos::RCP<Epetra_Vector> densityDiff = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    densityDiff->PutScalar(- restDensity_);
    densityDiff->Update(1.0, *densityn_, 1.0);

    // did we converge? Let's check the average
    densityDiff->Norm2(&densityErr);
    densityErr /= densityDiff->GlobalLength();
    if (densityErr <= correctDivergenceToll_)
    {
      break;
    }

    //std::cout << "The density error is = " << densityErr << std::endl;
    //std::cout << "Iteration " << ii << "/" << correctDensityIter_ << std::endl;
    //std::cin.get();

    // update the local ParticleMeshFreeData of the processor
    interHandler_->SetStateVector(densityn_, Density);

    // create and compute the pressure field required to have constant density
    Teuchos::RCP<Epetra_Vector> constDensityPressureAcc = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
    interHandler_->Inter_pvp_constDensityPressureAcc(constDensityPressureAcc, dt);

    // update the velocity
    veln_->Update(dt, *constDensityPressureAcc, 1.0);
  }

  std::cout << "Density corrected with " << ii << "/" << correctDensityIter_ << " iterations\n";
  std::cout << "Final density error = " << densityErr << std::endl;
}


/*--------------------------------------------------------------------------*
 | set the time step size and check that it is safe            katta 01/17  |
 *--------------------------------------------------------------------------*/
double PARTICLE::TimIntDivFree::SetDt()
{
  // decide the time step size
  double dt = DRT::Problem::Instance()->ParticleParams().get<double>("TIMESTEP");

  double radiusMin = 0;
  double velocityMax = 0;
  radiusn_->MinValue(&radiusMin);
  veln_->NormInf(&velocityMax);

  // small check on the radii
  if (radiusMin<0)
  {
    dserror("Negative radius... something strange is happening");
  }

  if (velocityMax > 0)
  {
    const double dt_CFL = PARTICLE::Utils::CFLcondition(radiusMin, velocityMax);

    switch (timeStepType_)
    {
    case INPAR::PARTICLE::Manual :
    {
      if (dt_CFL < dt)
      {
        std::cout << "Warning! some collision can be neglected! Use the Auto_CFL instead";
      }
      break;
    }
    case INPAR::PARTICLE::Auto_CFL :
    {
      if (dt_CFL < dt)
      {
        dt = dt_CFL;
      }
      break;
    }
    }
  }

  dt_->SetStep(0,dt);

  std::cout << "\n\n time step size: " << dt << std::endl << std::endl << std::endl;

  return dt;
}





