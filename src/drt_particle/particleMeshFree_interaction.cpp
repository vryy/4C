/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.cpp

\brief Particle-MeshFree interaction handling

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

/*
 References:
 Antoci2007: Numerical simulation of fluid–structure interaction by SPH, doi: 10.1016/j.compstruc.2007.01.002.
 Monaghan2005: Smoothed particle hydrodynamics, doi:10.1088/0034-4885/68/8/R01.
 Espanol2003: Smoothed dissipative particle dynamics, DOI: 10.1103/PhysRevE.67.026705.
 Akinci2013: http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf
*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particleMeshFree_interaction.H"
#include "particleMeshFree_weightFunction.H"
#include "particleMeshFree_rendering.H"
#include "particle_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "particle_timint_kickdrift.H"
#include "particle_heatSource.H"
#include "particle_utils.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 | constructor for particle-MeshFree interaction           katta 10/16  |
 *----------------------------------------------------------------------*/

PARTICLE::ParticleMeshFreeInteractionHandler::ParticleMeshFreeInteractionHandler(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams,
  const double initDensity,
  const double restDensity,
  const double refdensfac) :
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  myrank_(discret->Comm().MyPID()),
  weightFunctionHandler_(Teuchos::null),
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(particledynparams,"WALL_INTERACTION_TYPE")),
  freeSurfaceType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE")),
  WF_DIM_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunctionDim>(particledynparams,"WEIGHT_FUNCTION_DIM")),
  initDensity_(initDensity),
  restDensity_(restDensity),
  refdensfac_(refdensfac),
  periodic_length_(-1.0),
  min_pvp_dist_(1.0e10),
  min_pvw_dist_(1.0e10),
  min_pressure_(1.0e10),
  max_pressure_(-1.0e10)
{
  // checks
  if (particle_algorithm_->ExtParticleMat() == NULL)
    dserror("extParticleMat_ is empty");
  // extract material parameters
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
  {
    INPAR::PARTICLE::ExtendedGhosting extendedGhosting = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ExtendedGhosting>(particledynparams,"EXTENDED_GHOSTING");
    if(extendedGhosting != INPAR::PARTICLE::BdryParticleGhosting and extendedGhosting != INPAR::PARTICLE::AddLayerGhosting)
      dserror("Extended ghosting is required if boundary particles are applied!");
  }

  #ifndef PARTICLE_ONLYLAPLACETERM
  //Attention: The trace of a tensor as required to derive equation (28) from equation (25) in Espanol2003 depends
  //on the spatial dimension --> SPH approximations of the Laplace operator typically differ for different dimensions!
  switch (WF_DIM_)
  {
    case INPAR::PARTICLE::WF_3D :
    {
      convectionCoeff_ = 5.0*(extParticleMat->bulkViscosity_+extParticleMat->dynamicViscosity_ / 3.0);
      break;
    }
    case INPAR::PARTICLE::WF_2D :
    {
      convectionCoeff_ = 4.0*(extParticleMat->bulkViscosity_+extParticleMat->dynamicViscosity_ / 3.0);
      break;
    }
    case INPAR::PARTICLE::WF_1D :
    {
      convectionCoeff_ = 3.0*(extParticleMat->bulkViscosity_+extParticleMat->dynamicViscosity_ / 3.0);
      break;
    }
    default :
    {
      dserror("Only the problem dimensions 1, 2 and 3 are possible!");
    }
  }

  diffusionCoeff_ = 5.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
  #else
  diffusionCoeff_ = 2.0 * extParticleMat->dynamicViscosity_;
  convectionCoeff_= 0.0;
  if(extParticleMat->bulkViscosity_>0)
    dserror("Bulk Viscosity not considered here!");
  #endif

  // checks
  if (diffusionCoeff_<0)
  {
    dserror("The diffusion coefficient is negative! The following equation should hold: 5*dynamicViscosity >= 3*bulkViscosity");
  }

  if (convectionCoeff_<0)
  {
    dserror("The convection coefficient is negative! Are you sure that the dynamic viscosity and the bulk modulus are positive?");
  }

  artificialViscosity_ = extParticleMat->bulkViscosity_+extParticleMat->artificialViscosity_;

  surfaceTension_ = extParticleMat->surfaceTension_;
  staticContactAngle_ = extParticleMat->staticContactAngle_;

  // set the correct WeightFunction
  switch (DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(particledynparams,"WEIGHT_FUNCTION"))
  {
    case INPAR::PARTICLE::CubicBspline :
    {
      weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_CubicBspline(WF_DIM_));
      break;
    }
    case INPAR::PARTICLE::QuinticBspline :
    {
      weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_QuinticBspline(WF_DIM_));
      break;
    }
    case INPAR::PARTICLE::SqrtHyperbola :
    {
      weightFunctionHandler_ = Teuchos::rcp(new WeightFunction_SqrtHyperbola(WF_DIM_));
      break;
    }
    case INPAR::PARTICLE::HyperbolaNoRsz :
    {
      weightFunctionHandler_ = Teuchos::rcp(new WeightFunction_HyperbolaNoRsz);
      break;
    }
  }

  if(DRT::INPUT::IntegralValue<int>(particledynparams,"SOLVE_THERMAL_PROBLEM")==false and extParticleMat->initTemperature_<extParticleMat->transitionTemperature_)
    dserror("Pure mechanical problems, i.e. SOLVE_THERMAL_PROBLEM==No, are currently only considered as pure fluid problems. "
        "To remain consistent, the initial temperature has to be higher than the transition temperature!");

  //Check if periodic boundary conditions are applied
  BuildPeriodicBC();

  // initialize rendering handler
  const INPAR::PARTICLE::RenderingType renderingType = DRT::INPUT::IntegralValue<INPAR::PARTICLE::RenderingType>(DRT::Problem::Instance()->ParticleParams(),"RENDERING");
  if (renderingType == INPAR::PARTICLE::StandardRendering or renderingType == INPAR::PARTICLE::NormalizedRendering)
  {
    Teuchos::RCP<Rendering> rendering = Teuchos::rcp(new PARTICLE::Rendering(particle_algorithm_, Teuchos::rcp(this,false), weightFunctionHandler_));
    particle_algorithm_->SetRendering(rendering);

    if (particle_algorithm_->GetRendering() == Teuchos::null)
      dserror("Rendering handler not set correctly in particle algorithm!");
  }

#ifndef PARTICLE_BOUNDARYDENSITY
  if(PARTICLE_REINITSHIFT!=2)
    dserror("If PARTICLE_REINITSHIFT != 2, the flag PARTICLE_BOUNDARYDENSITY is required!");
#endif
}

/*----------------------------------------------------------------------*
 | build periodic boundary conditions                       meier 03/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::BuildPeriodicBC()
{

  std::istringstream periodicbc(Teuchos::getNumericStringParameter(
      DRT::Problem::Instance()->BinningStrategyParams(),"PERIODICONOFF"));

  bool periodic_bcs=false;

  // loop over all spatial directions
  for(int dim=0; dim<3; ++dim)
  {
    int val = -1;
    if (periodicbc >> val)
    {
      if( val )
      {
        if(dim==0)
        {
          if(myrank_ == 0)
            std::cout << "INFO: Periodic boundary conditions in x-direction are applied!" << std::endl;

          periodic_bcs=true;
        }
        else
          dserror("Periodic boundary conditions only possible in x-direction so far!");
      }
    }
    else
    {
      dserror("Enter three values to specify each direction as periodic or non periodic (y and z value has to be zero). Fix input file ...");
    }
  }

  if(periodic_bcs)
  {
    LINALG::Matrix<3,2> box(true);

    // get bounding box specified in the input file
    box.PutScalar(1.0e12);
    std::istringstream xaabbstream( Teuchos::getNumericStringParameter(
        DRT::Problem::Instance()->BinningStrategyParams(),"BOUNDINGBOX") );
    for( int col = 0; col < 2; ++col )
    {
      for( int row = 0; row < 3; ++row )
      {
        double value = 1.0e12;
        if( xaabbstream >> value )
          box( row, col ) = value;
        else
          dserror(" Specify six values for bounding box in three dimensional problem."
                  " Fix your input file.");
      }
    }

    if(fabs(box(0,0)+box(0,1))>0 or box(0,0)>box(0,1))
      dserror("Left and Right periodic boundary have to be of the form x_b=-periodic_length_/2 and x_b=periodic_length_/2!");
    else
      periodic_length_=2.0*fabs(box(0,1));
  }

  return;
}


/*----------------------------------------------------------------------*
 | set up internal variables for future computations       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Init(
    const int step,
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn)
{
  // check
  if (colParticles_.size() != 0 or colFADParticles_.size() != 0)
  {
    dserror("you did not call Clear before Init, colParticles_ is not empty");
  }

  // set up the local data storage and fill it with the state vectors
  if (!neighbours_p_.empty() || !neighbours_hs_.empty() || !overlappingneighbours_p_.empty())
  {
    std::cout << "The neighbours have memory (They are not empty)!\n";
    std::cout << "However, lid row/col ids were not updated (because not yet updated). It is safe only when not parallel\n";
    std::cin.get();
  }
  InitColParticles();

  // set up positions and radii to set up the neighbours
  if(disn!=Teuchos::null)
    SetStateVector(disn, PARTICLE::Dis);

  if(radiusn!=Teuchos::null)
    SetStateVector(radiusn, PARTICLE::Radius);

  // keep going with the remaining state vectors
  if(veln!=Teuchos::null)
    SetStateVector(veln, PARTICLE::Vel);
  if(mass!=Teuchos::null)
    SetStateVector(mass, PARTICLE::Mass);
  if(specEnthalpyn!=Teuchos::null)
    SetStateVector(specEnthalpyn, PARTICLE::SpecEnthalpy);

  // set up the neighbours
  AddNewNeighbours(step);
}

/*----------------------------------------------------------------------*
 | set up internal variables for future computations       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Init(
    const int step,
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> temperature)
{
  Init(step, disn, veln, radiusn, mass, specEnthalpyn);

  // set the other state vectors
  SetStateVector(temperature, PARTICLE::Temperature);
}

/*----------------------------------------------------------------------*
 | set up internal variables for future computations       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Init(
    const int step,
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<const Epetra_Vector> densityn,
    Teuchos::RCP<const Epetra_Vector> pressure)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::Init");

  Init(step, disn, veln, radiusn, mass, specEnthalpyn, temperature);

  // set the other state vectors
  SetStateVector(densityn, PARTICLE::Density);
  SetStateVector(pressure, PARTICLE::Pressure);
}


/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitColParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::InitColParticles");
  const int numcolelements = discret_->NodeColMap()->NumMyElements();

  colParticles_.resize(numcolelements);
  colFADParticles_.resize(numcolelements);
  boundaryparticles_.clear();

  for (int lidNodeCol=0; lidNodeCol<numcolelements; ++lidNodeCol)
  {
    DRT::Node *particle = discret_->lColNode(lidNodeCol);
    std::vector<int> lm;
    lm.reserve(3);
    discret_->Dof(particle, lm);

    bool boundaryparticle = false;
    if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    {
      PARTICLE::ParticleNode* particleNode = dynamic_cast<PARTICLE::ParticleNode*>(particle);
      if (particleNode == NULL)
        dserror("Dynamic cast to ParticleNode failed");
      boundaryparticle = particleNode->Is_bdry_particle();
    }

    colParticles_[lidNodeCol] = ParticleMF(particle->Id(), particle->Owner(), lm, boundaryparticle);
    colFADParticles_[lidNodeCol] = Teuchos::rcp(new ParticleFAD(particle->Id()));

    if(boundaryparticle)
    {
      boundaryparticles_[particle->Id()]=&colParticles_[lidNodeCol];
    }
  }
}

/*-------------------------------------------------------------------------------------------------*
 | accelerations, modified pressures, velocities and energy for boundary particles    meier 02/17  |
 *-------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitBoundaryData(
        Teuchos::RCP<const Epetra_Vector> accn,
        const LINALG::Matrix<3,1>& g,
        double &bpintergy)
{
  //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    dserror("Boundary particles are only considered in the context of pure fluid problems so far!");

  SetStateVector(accn, PARTICLE::Acc);

  //Clear energy in the beginning
  bpintergy=0.0;

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;

    ParticleMF *particle_i = boundaryparticles_[id_i];

    LINALG::Matrix<3,1> sumjVjWij(true);
    LINALG::Matrix<3,1> sumjVDiffjWij(true);
    double sumjWij(0.0);

    double sumjpjWij(0.0);
    LINALG::Matrix<3,1> sumjrhojrijWij(true);

    LINALG::Matrix<3,1> gminusacci(true);
    gminusacci.Update(1.0,g,1.0);
    gminusacci.Update(-1.0,particle_i->boundarydata_.acc_,1.0);

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {

      int id_j = jj->first;
      int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      //sum up the numerator of Ref. Adami2012, Eq.(22)
      sumjVjWij.Update(interData_ij.w_ij_,particle_j.vel_,1.0);
      //sum up the denominator of Ref. Adami2012, Eq.(22)
      sumjWij+=interData_ij.w_ij_;
      //sum up the left sum of numerator of Ref. Adami2012, Eq.(27)
      sumjpjWij+=particle_j.pressure_*interData_ij.w_ij_;
      //sum up the right sum of numerator of Ref. Adami2012, Eq.(27)
      sumjrhojrijWij.Update(interData_ij.rRelNorm2_*particle_j.density_*interData_ij.w_ij_,interData_ij.rRelVersor_ij_,1.0);

      //sum up difference velocity required for boundary treatment of tensor A in transport velocity formulation according to
      //Adami2013. Note: In the original formulation Adami2013, no boundary treatment of tensor A has been applied.
      sumjVDiffjWij.Update(interData_ij.w_ij_,particle_j.velConv_,1.0);
      sumjVDiffjWij.Update(-interData_ij.w_ij_,particle_j.vel_,1.0);
    }

    //Only set valuces in case there are interacting (=close enough) fluid particles
    if(sumjWij>0)
    {
      //Update boundary particle velocity required for viscous forces, see Ref. Adami2012, Eq.(23)
      LINALG::Matrix<3,1> tildeVi(true);
      tildeVi.Update(1.0/sumjWij,sumjVjWij,0.0);

      LINALG::Matrix<3,1> vw(true);
      vw.Update(2.0,particle_i->vel_,1.0);
      vw.Update(-1.0,tildeVi,1.0);

      //In case of a no-slip boundary condition, we apply a quasi-linear extrapolation according to Ref. Adami2012, Eq.(23)
      particle_i->boundarydata_.velModVisc_.Update(1.0,vw,0.0);

      //Update boundary particle pressure required for pressure forces, see Ref. Adami2012, Eq.(23)
      double dotProduct = gminusacci.Dot(sumjrhojrijWij);
      double pi=(sumjpjWij+dotProduct)/sumjWij;
      particle_i->boundarydata_.pressureMod_=pi;

      //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
      double density_i=PARTICLE::Utils::Pressure2Density(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),
          restDensity_,refdensfac_,pi,particle_algorithm_->ExtParticleMat()->exponent_);
      particle_i->boundarydata_.densityMod_=density_i;

      //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
      bpintergy+=PARTICLE::Utils::Density2Energy(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),
          density_i, restDensity_, refdensfac_, particle_i->mass_);
    }
    else
    {
      particle_i->boundarydata_.pressureMod_=-1000; //-1000 means there are no relevant fluid particle neighbors for this boundary particle
      particle_i->boundarydata_.densityMod_=-1000; //-1000 means there are no relevant fluid particle neighbors for this boundary particle
    }
  }
}


/*----------------------------------------------------------------------*
 | clear data, keep memory                                 katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Clear()
{
  // erase colParticles_. keep the memory
  colParticles_.clear();
  colFADParticles_.clear();
  // erase neighbours keep memory
  neighbours_p_.clear();
  neighbours_hs_.clear();
  overlappingneighbours_p_.clear();
  neighbours_bp_.clear();
  neighbours_fsp_.clear();
  boundaryparticles_.clear();
  freesurfaceparticles_.clear();
}


/*--------------------------------------------------------------------------------------*
 | clear data, keep the memory and keep interactions up to step-memory     katta 12/16  |
 *--------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Clear(const int step, const int memory)
{
  // erase colParticles_. keep the memory
  colParticles_.clear();
  colFADParticles_.clear();
  // erase neighbours keep memory
  neighbours_hs_.clear();
  neighbours_bp_.clear();
  neighbours_fsp_.clear();
  boundaryparticles_.clear();
  freesurfaceparticles_.clear();

  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      if (jj->second.step_ + memory <= step)
      {
        jj = neighbours_p_[lidNodeRow_i].erase(jj);
      }
    }
  }

  #ifdef PARTICLE_OVERLAPPINGNEIGHBORS
  for (unsigned int lidNodeCol_i = 0; lidNodeCol_i != overlappingneighbours_p_.size(); ++lidNodeCol_i)
  {
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = overlappingneighbours_p_[lidNodeCol_i].begin(); jj != overlappingneighbours_p_[lidNodeCol_i].end(); ++jj)
    {
      if (jj->second.step_ + memory <= step)
      {
        jj = overlappingneighbours_p_[lidNodeCol_i].erase(jj);
      }
    }
  }
  #endif
}

/*----------------------------------------------------------------------*
 | set colVectors in the local data structs                katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVector(Teuchos::RCP<const Epetra_Vector> stateVector, const PARTICLE::StateVectorType svt)
{
  // checks
  if (stateVector == Teuchos::null)
  {
    dserror("the state vector is empty");
  }

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVector");

  /// miraculous transformation into column vector... ///
  Teuchos::RCP<Epetra_Vector> stateVectorCol;
  switch (svt)
  {
    // dof based vectors
    case PARTICLE::Dis :
    case PARTICLE::Vel :
    case PARTICLE::VelConv :
    case PARTICLE::Acc :
    {
      stateVectorCol = LINALG::CreateVector(*discret_->DofColMap(),false);
      break;
    }
    // node based vectors
    default :
    {
      stateVectorCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
      break;
    }
  }
  LINALG::Export(*stateVector ,*stateVectorCol);

  // fill particleData_
  for (int lidNodeCol=0; lidNodeCol<discret_->NodeColMap()->NumMyElements(); ++lidNodeCol)
  {
    ParticleMF& data = colParticles_[lidNodeCol];

    switch (svt)
    {
      // dof based vectors
      case PARTICLE::Dis :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.dis_, data.lm_);

        break;
      }
      case PARTICLE::Vel :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.vel_, data.lm_);
        break;
      }
      case PARTICLE::VelConv :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.velConv_, data.lm_);
        break;
      }
      case PARTICLE::Acc :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.boundarydata_.acc_, data.lm_);
        break;
      }
      // node based vectors
      case PARTICLE::Radius :
      {
        data.radius_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::Density :
      {
        data.density_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::DensityDot :
      {
        data.densityDot_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::SpecEnthalpy :
      {
        data.specEnthalpy_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::Mass :
      {
        data.mass_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::Pressure :
      {
        data.pressure_ = (*stateVectorCol)[lidNodeCol];
        if(data.boundarydata_.boundaryparticle_==false and data.pressure_ > max_pressure_)
          max_pressure_=data.pressure_;
        if(data.boundarydata_.boundaryparticle_==false and data.pressure_ < min_pressure_)
          min_pressure_=data.pressure_;
        break;
      }
      case PARTICLE::Temperature :
      {
        data.temperature_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::ColorField :
      {
        data.freesurfacedata_.color_field_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::DensitySum :
      {
        data.freesurfacedata_.density_sum_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      default :
      {
        dserror("StateVector type unknown or wrong stateVector format inserted");
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours(const int step)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours");

  // resize the vectors
  const int numRowParticles = discret_->NodeRowMap()->NumMyElements();
  neighbours_p_.resize(numRowParticles);
  neighbours_hs_.resize(numRowParticles);

  #ifdef PARTICLE_OVERLAPPINGNEIGHBORS
  const int numColParticles = discret_->NodeColMap()->NumMyElements();
  overlappingneighbours_p_.resize(numColParticles);
  #endif

  // bin checker
  std::set<int> examinedbins;

  //***********loop over the (real/fluid) particles (no superpositions)****************************************************************************************************************
  const int numrowparticles = discret_->NodeRowMap()->NumMyElements();
  for(int rowPar_i=0; rowPar_i<numrowparticles; ++rowPar_i)
  {
    // extract the particle
    DRT::Node *currparticle = discret_->lRowNode(rowPar_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // list of particles in Bin
    std::list<DRT::Node*> neighboursLinf_p;

    // list of walls that border on the CurrentBin
    boost::unordered_map<int, DRT::Element*> neighboursLinf_w;

    // list of heat sources that border on the CurrentBin
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::rcp(new boost::unordered_map<int , Teuchos::RCP<HeatSource> >);

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p, &neighboursLinf_w, neighboursLinf_hs);

    // do some checks
    if(neighboursLinf_hs->size()!=0)
    {
      if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip)
        dserror("Combination of boundary particles and heat sources not possible so far!");

      if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM")==false)
        dserror("Heat sources are defined but thermal problem is not solved. Is that your intention?");
    }

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      //std::cout << "particle_i.gid_: " << particle_i.gid_ << std::endl;

      AddNewNeighbours_p(particle_i, neighboursLinf_p, step);

      if(neighboursLinf_hs->size()!=0)
        AddNewNeighbours_hs(particle_i, neighboursLinf_hs);
    }
  }

  //*************loop over the boundary particles (for boundary particles we need superposition --> column map)**************************************************
  // clear bin checker
  examinedbins.clear();

  // loop over the boundary particles
  for (std::map<int,ParticleMF* >::iterator ii = boundaryparticles_.begin(); ii!=boundaryparticles_.end(); ++ii)
  {
    int id_i = ii->first;
    // extract the node underlying the particle
    DRT::Node *currparticle = discret_->gNode(id_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins2_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // list of particles in Bin
    std::list<DRT::Node*> neighboursLinf_p;

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(particle_i.boundarydata_.boundaryparticle_==true)
      {
        AddNewNeighbours_bp(particle_i, neighboursLinf_p, step);
      }
    }
  }


  //***********loop over the (real/fluid) particles (overlapping vector)****************************************************************************************************************
  #ifdef PARTICLE_OVERLAPPINGNEIGHBORS// clear bin checker
  examinedbins.clear();

  const int numcolparticles = discret_->NodeColMap()->NumMyElements();

  for(int colPar_i=0; colPar_i<numcolparticles; ++colPar_i)
  {
    // extract the particle
    DRT::Node *currparticle = discret_->lColNode(colPar_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // list of particles in Bin
    std::list<DRT::Node*> neighboursLinf_p;

    // list of walls that border on the CurrentBin
    boost::unordered_map<int, DRT::Element*> neighboursLinf_w;

    // list of heat sources that border on the CurrentBin
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::rcp(new boost::unordered_map<int , Teuchos::RCP<HeatSource> >);

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p);



    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {

      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];


      AddNewNeighbours_op(particle_i, neighboursLinf_p, step);
    }
  }
  #endif
}


/*----------------------------------------------------------------------*
 | set the neighbours - particles                          katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_p(
    const ParticleMF& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p,
    const int step)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_p");

  // self-neighbours not allowed
  // insert the interaction only if meaningful

  // loop over the neighbours particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

    //Attention: In the current framework, two different neighbor vectors neighbours_p_ and overlappingneighbours_p_ are filled,
    //neighbours_p_ (node row map) takes advantage of the symmetry of contact interaction and only stores the pair ij but not the pair ji for i<j
    //For debugging purposes and non-symmetric interaction laws additionally the fully overlapping vector vectorneighbours_p_ is filled.
    //The latter contains both pairs ij and ji and is furthermore based on a node column map
    //(which will be relevant for the calculation of linearizations where "the neighbors of neighbors" will be required).
    //Moreover, the vector overlappingneighbours_p_ also contains auto/self interaction pairs ii!
    if(particle_j.owner_ != myrank_ || (particle_i.gid_ < particle_j.gid_))
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      LINALG::Matrix<3,1> dis_i=particle_i.dis_;
      LINALG::Matrix<3,1> dis_j=particle_j.dis_;

      //In case both particles i and j are close to two opposite periodic boundaries, the position of particle j is shifted correspondingly.
      //It is sufficient to only shift dis_j since the terms evaluated in the following do only depend on the relative distance between two
      //particles but not on the absolute positions.
      ShiftPeriodicBoundaryPair(dis_i,dis_j,particle_i.radius_,particle_j.radius_);

      rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

      //Only evaluate pairs that are close enough.
      //TODO: In future applications (e.g. implicit time integration), where neighbors are not searched for every new configuration
      //this criterion might be to strict!!!
      if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
        continue;

      #ifdef PARTICLE_TENSILESAFETYFAC
        //Check that the particle distance does not become to small (danger of tensile instabilities).
        //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
        double initAverageDist = PARTICLE::Utils::Volume2EffDist(std::max(particle_i.mass_,particle_j.mass_)/initDensity_,WF_DIM_);
        if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
          dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
      #endif

      //check for minimal distance
      if(rRelNorm2<min_pvp_dist_)
        min_pvp_dist_=rRelNorm2;

      const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
      const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
      const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

      double w_ji = 0;
      double dw_ji = 0;
      double ddw_ji = 0;

      //if (particle_j.owner_ == myrank_)
      //{
        if (particle_i.radius_ == particle_j.radius_)
        {
            w_ji = w_ij;
           dw_ji = dw_ij;
          ddw_ji = ddw_ij;
        }
        else
        {
            w_ji = weightFunctionHandler_->W(rRelNorm2, particle_j.radius_);
           dw_ji = weightFunctionHandler_->DW(rRelNorm2, particle_j.radius_);
          ddw_ji = weightFunctionHandler_->DDW(rRelNorm2, particle_j.radius_);
        }
      //}

      // push_back
      LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
      rRelVersor_ij.Scale(1/rRelNorm2);
      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
      (neighbours_p_[lidNodeRow_i])[lidNodeCol_j] = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          step,
          w_ij,
          w_ji,
          dw_ij,
          dw_ji,
          ddw_ij,
          ddw_ji);
    }
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - particles (overlapping vector)     meier 03/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_op(
    const ParticleMF& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p,
    const int step)
{
  // self-neighbours not allowed
  // insert the interaction only if meaningful

  // loop over the neighbours particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {

    const int lidNodeCol_j = (*jj)->LID();
    const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

    //Attention: In the current framework, two different neighbor vectors neighbours_p_ and overlappingneighbours_p_ are filled,
    //neighbours_p_ (node row map) takes advantage of the symmetry of contact interaction and only stores the pair ij but not the pair ji for i<j
    //For debugging purposes and non-symmetric interaction laws additionally the fully overlapping vector vectorneighbours_p_ is filled.
    //The latter contains both pairs ij and ji and is furthermore based on a node column map
    //(which will be relevant for the calculation of linearizations where "the neighbors of neighbors" will be required).
    //Moreover, the vector overlappingneighbours_p_ also contains auto/self interaction pairs ii!

    double rRelNorm2 = 0.0;
    LINALG::Matrix<3,1> rRelVersor_ij(true);
    double w_ij = 0.0;
    double dw_ij = 0.0;
    double ddw_ij = 0.0;

    const int lidNodeCol_i = discret_->NodeColMap()->LID(particle_i.gid_);
    if(lidNodeCol_i!=lidNodeCol_j)//standard case: interaction ij with i != j
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      LINALG::Matrix<3,1> dis_i=particle_i.dis_;
      LINALG::Matrix<3,1> dis_j=particle_j.dis_;

      //In case both particles i and j are close to two opposite periodic boundaries, the position of particle j is shifted correspondingly.
      //It is sufficient to only shift dis_j since the terms evaluated in the following do only depend on the relative distance between two
      //particles but not on the absolut positions.
      ShiftPeriodicBoundaryPair(dis_i,dis_j,particle_i.radius_,particle_j.radius_);

      rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
      rRelNorm2 = rRel_ij.Norm2();

      //Only evaluate pairs that are close enough.
      //TODO: In future applications (e.g. implicit time integration), where neighbors are not searched for every new configuration
      //this criterion might be to strict!!!
      if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
        continue;

      rRelVersor_ij.Update(1.0,rRel_ij,0.0);
      rRelVersor_ij.Scale(1/rRelNorm2);

      w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
      dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
      ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

      #ifdef PARTICLE_TENSILESAFETYFAC
        //Check that the particle distance does not become to small (danger of tensile instabilities).
        //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
        double initAverageDist = PARTICLE::Utils::Volume2EffDist(std::max(particle_i.mass_,particle_j.mass_)/initDensity_,WF_DIM_);
        if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
          dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
      #endif
    }
    else //auto/self interaction ii
    {
      rRelNorm2 = 0.0;
      rRelVersor_ij.Clear();
      w_ij = weightFunctionHandler_->W0(particle_i.radius_);
      dw_ij = 0.0;
      ddw_ij = weightFunctionHandler_->DDW0(particle_i.radius_);
    }

    //For the overlapping vector, only paris ij but not ji are required!
    double w_ji = -1000;
    double dw_ji = -1000;
    double ddw_ji = -1000;

    (overlappingneighbours_p_[lidNodeCol_i])[lidNodeCol_j] = InterDataPvP(
        rRelVersor_ij,
        rRelNorm2,
        step,
        w_ij,
        w_ji,
        dw_ij,
        dw_ji,
        ddw_ij,
        ddw_ji);
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - boundary particles                 meier 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_bp(
    const ParticleMF& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p,
    const int step)
{
  // self-neighbours not allowed, boundary particles not allowed as neighbors

  // loop over the neighbours_ particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

    //boundary particles not allowed as neighbors
    if(particle_j.boundarydata_.boundaryparticle_==false)
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      LINALG::Matrix<3,1> dis_i=particle_i.dis_;
      LINALG::Matrix<3,1> dis_j=particle_j.dis_;

      //In case both particles i and j are close to two opposite periodic boundaries, the position of particle j is shifted correspondingly.
      //It is sufficient to only shift dis_j since the terms evaluated in the following do only depend on the relative distance between two
      //particles but not on the absolut positions.
      ShiftPeriodicBoundaryPair(dis_i,dis_j,particle_i.radius_,particle_j.radius_);

      rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

      //Only evaluate pairs that are close enough.
      //TODO: In future applications (e.g. implicit time integration), where neighbors are not searched for every new configuration
      //this criterion might be to strict!!!
      if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
        continue;

      #ifdef PARTICLE_TENSILESAFETYFAC
        //Check that the particle distance does not become to small (danger of tensile instabilities).
        //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
        double initAverageDist = PARTICLE::Utils::Volume2EffDist(std::max(particle_i.mass_,particle_j.mass_)/initDensity_,WF_DIM_);
        if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
          dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
      #endif

      //check for minimal distance
       if(rRelNorm2<min_pvw_dist_)
         min_pvw_dist_=rRelNorm2;

      const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
      const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
      const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

      double w_ji = 0;
      double dw_ji = 0;
      double ddw_ji = 0;

      if (particle_i.radius_ == particle_j.radius_)
      {
          w_ji = w_ij;
         dw_ji = dw_ij;
        ddw_ji = ddw_ij;
      }
      else
      {
          w_ji = weightFunctionHandler_->W(rRelNorm2, particle_j.radius_);
         dw_ji = weightFunctionHandler_->DW(rRelNorm2, particle_j.radius_);
        ddw_ji = weightFunctionHandler_->DDW(rRelNorm2, particle_j.radius_);
      }

    // push_back
    LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
    rRelVersor_ij.Scale(1/rRelNorm2);
    (neighbours_bp_[particle_i.gid_])[particle_j.gid_] = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          step,
          w_ij,
          w_ji,
          dw_ij,
          dw_ji,
          ddw_ij,
          ddw_ji);
    }
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - free-surface particles             meier 08/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_fsp(
    const ParticleMF& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p,
    const bool mark_fs_indirect,
    const int step)
{
  //TODO: We do not really need the variable step so far. Thus, it is set to 1 by default so far.
  // self-neighbours not allowed, boundary particles not allowed as neighbors

  // loop over the neighbours_ particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    ParticleMF& particle_j = colParticles_[lidNodeCol_j];

    if(particle_i.gid_==particle_j.gid_)
      continue;

    // create the data that we have to insert
    LINALG::Matrix<3,1> rRel_ij;
    LINALG::Matrix<3,1> dis_i=particle_i.dis_;
    LINALG::Matrix<3,1> dis_j=particle_j.dis_;

    if(periodic_length_>0)
      dserror("Periodic free surface flows not possible so far!");

    rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
    const double rRelNorm2 = rRel_ij.Norm2();

    //Only evaluate pairs that are close enough.
    //TODO: In future applications (e.g. implicit time integration), where neighbors are not searched for every new configuration
    //this criterion might be to strict!!!
    if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
      continue;

    #ifdef PARTICLE_TENSILESAFETYFAC
      //Check that the particle distance does not become to small (danger of tensile instabilities).
      //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
      double initAverageDist = PARTICLE::Utils::Volume2EffDist(std::max(particle_i.mass_,particle_j.mass_)/initDensity_,WF_DIM_);
      if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
        dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
    #endif

    const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
    const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
    const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

    double w_ji = 0;
    double dw_ji = 0;
    double ddw_ji = 0;

    if (particle_i.radius_ == particle_j.radius_)
    {
        w_ji = w_ij;
       dw_ji = dw_ij;
      ddw_ji = ddw_ij;
    }
    else
    {
        w_ji = weightFunctionHandler_->W(rRelNorm2, particle_j.radius_);
       dw_ji = weightFunctionHandler_->DW(rRelNorm2, particle_j.radius_);
      ddw_ji = weightFunctionHandler_->DDW(rRelNorm2, particle_j.radius_);
    }

  // push_back
  LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
  rRelVersor_ij.Scale(1/rRelNorm2);
  (neighbours_fsp_[particle_i.gid_])[particle_j.gid_] = InterDataPvP(
        rRelVersor_ij,
        rRelNorm2,
        step,
        w_ij,
        w_ji,
        dw_ij,
        dw_ji,
        ddw_ij,
        ddw_ji);

  //mark particle_j as indirect free surface particle in case in has not yet been marked as direct free-surface particle
  if(mark_fs_indirect and particle_j.freesurfacedata_.freesurfaceparticletype_!=FS_DIRECT and particle_j.boundarydata_.boundaryparticle_==false)
    particle_j.freesurfacedata_.freesurfaceparticletype_=FS_INDIRECT;
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - heat sources                       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_hs(
    const ParticleMF& particle_i,
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > >& neighboursLinf_hs)
{
  boost::unordered_map<int , Teuchos::RCP<HeatSource> >::const_iterator ii;
  for(ii = neighboursLinf_hs->begin(); ii != neighboursLinf_hs->end();  ++ii)
  {
    const Teuchos::RCP<HeatSource> hs = ii->second;
    if (hs->minVerZone_[0]<=particle_i.dis_(0) &&
        hs->minVerZone_[1]<=particle_i.dis_(1) &&
        hs->minVerZone_[2]<=particle_i.dis_(2) &&
        hs->maxVerZone_[0]>=particle_i.dis_(0) &&
        hs->maxVerZone_[1]>=particle_i.dis_(1) &&
        hs->maxVerZone_[2]>=particle_i.dis_(2))
    {
      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
      (neighbours_hs_[lidNodeRow_i])[hs->id_] = hs;
    }
  }
}

/*---------------------------------------------------------------------------------*
 | Apply coordinate shift in case of periodic boundary conditions     meier 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::ShiftPeriodicBoundaryPair(
    const LINALG::Matrix<3,1>& dis_1,
    LINALG::Matrix<3,1>& dis_2,
    const double& radius_1,
    const double& radius_2)
{
  // no periodic boundary conditions
  if( periodic_length_ <= 0 )
    return;

  //We dont want to have particles that can interact directly AND across the periodic boundary
  //(which would be possible in case of periodic_length_<4*particle_i.radius_)
  //Update 05/24/2017: Actually, the factor 4 could be replaced by 2 in the SPH case since interaction
  //is only possible for distances < radius (and not for distance < 2 radius as in DEM)!
  if(periodic_length_<4.0*std::max(radius_1,radius_2))
    dserror("The length periodic_length_ of the periodic box has to be larger than four times the particle radius!");

  double x_coord_1=dis_1(0);
  double x_coord_2=dis_2(0);

  if(x_coord_1>periodic_length_/2 or x_coord_1<-periodic_length_/2 or x_coord_2>periodic_length_/2 or x_coord_2<-periodic_length_/2)
    dserror("Particle outside the periodic domain!");

  //Within this method, only the position of particle 2 (required for interaction with particle 1) is shifted, while particle 1 remains untouched.
  if(x_coord_1>=periodic_length_/2-std::max(radius_1,radius_2) and x_coord_2<=-periodic_length_/2+std::max(radius_1,radius_2))
    dis_2(0)+=periodic_length_;
  else if(x_coord_2>=periodic_length_/2-std::max(radius_1,radius_2) and x_coord_1<=-periodic_length_/2+std::max(radius_1,radius_2))
    dis_2(0)-=periodic_length_;

  return;
}

/*----------------------------------------------------------------------*
 | evaluate density - pvp                                  katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_densityDot(
    const Teuchos::RCP<Epetra_Vector> densityDotn,
    const double extMulti)
{
  //checks
  if (densityDotn == Teuchos::null)
  {
    dserror("densityDotn is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    double sumj_densityDotn_ij=0.0;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //Boundary particles are not considered here
      if(particle_j.boundarydata_.boundaryparticle_==true)
        continue;

      // --- extract and compute general data --- //

      LINALG::Matrix<3,1> vRel_ij;

      bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");

      if(transport_velocity==false)
        vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);
      else
        vRel_ij.Update(1.0, particle_i.velConv_, -1.0, particle_j.velConv_);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // assemble and write
        LINALG::Matrix<3,1> gradW = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ij_);
        double densityDotn_ij = gradW.Dot(vRel_ij) * particle_j.mass_;
        densityDotn_ij *= extMulti;
        sumj_densityDotn_ij+=densityDotn_ij;
      }

      // write on particle j if appropriate specializing the quantities
      if ( interData_ij.dw_ji_ != 0)
      {
        // assemble and write
        // actio = - reaction twice... no change in sign required (vRel and rRel both change sign)
        LINALG::Matrix<3,1> gradW = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ji_);
        double densityDotn_ji = gradW.Dot(vRel_ij) * particle_i.mass_;
        densityDotn_ji *= extMulti;
        LINALG::Assemble(*densityDotn, densityDotn_ji, particle_j.gid_, particle_j.owner_);
      }
    }//sum over j
    //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*densityDotn, sumj_densityDotn_ij, particle_i.gid_, particle_i.owner_);
  }//sum over i

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    Inter_bpvp_densityDot(densityDotn);
}


/*----------------------------------------------------------------------*
 | evaluate density - bpvp                                 meier 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_bpvp_densityDot(
    const Teuchos::RCP<Epetra_Vector> densityDotn,
    const double extMulti)
{
  //checks
  if (densityDotn == Teuchos::null)
  {
    dserror("densityDotn is empty");
  }

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;
    ParticleMF *particle_i = boundaryparticles_[id_i];

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {
      int id_j = jj->first;
      int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      //Atenttion:  particle_i is allways a boundary particle and particle_j is allways a fluid particle.
      //            Forces are only assembled to the fluid particles particle_j!
      //            Everything else is very similar to Inter_pvp_densityDot().

      // --- extract and compute general data --- //

      LINALG::Matrix<3,1> vRel_ij;
      bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");
      if(transport_velocity==false)
        vRel_ij.Update(1.0, particle_i->vel_, -1.0, particle_j.vel_);
      else
        vRel_ij.Update(1.0, particle_i->velConv_, -1.0, particle_j.velConv_);

      // write on particle j if appropriate specializing the quantities
      if ( interData_ij.dw_ji_ != 0)
      {
        // assemble and write
        // actio = - reaction twice... no change in sign required (vRel and rRel both change sign)
        LINALG::Matrix<3,1> gradW = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ji_);
        double densityDotn_ji = gradW.Dot(vRel_ij) * particle_i->mass_;
        densityDotn_ji *= extMulti;
        LINALG::Assemble(*densityDotn, densityDotn_ji, particle_j.gid_, particle_j.owner_);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | evaluate acceleration - pvp                             katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_acc_var1(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> trvl_acc,
    const Teuchos::RCP<Epetra_Vector> acc_A,
    const double time,
    const double extMulti)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_acc");

  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }

  // Modified particle transport velocities either due to constant background pressure according to Adami2013 or due to XSPH according to Monaghan1989
  // The terms associated with the constant background pressure are similar to the real pressure terms
  // but with the pressures p_i and p_j replaced by the background pressure. All the related variables
  // and terms are marked by a prefix p0, e.g. p0_momentum_ij etc.. All terms related to XSPH are marked by a prefix xsph
  double background_pressure=-1.0;
  double xsph_dampfac=-1.0;
  double xsph_stiffac=-1.0;
  bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");
  if (trvl_acc != Teuchos::null)
  {
    background_pressure=DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE");
    xsph_dampfac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_DAMPFAC");
    xsph_stiffac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_STIFFAC");

    if(background_pressure<0.0 and xsph_dampfac<0.0 and xsph_stiffac<0.0)
      dserror("Negative values of background pressure, xsph_dampfac and xsph_stifffac! Why is the vector accmod handed in for this case? Set flag TRANSPORT_VELOCITY to no!");

    if(transport_velocity==false)
      dserror("Vector accmod handed in even though TRANSPORT_VELOCITY is set to no?!?");

    if((xsph_dampfac<0.0 and xsph_stiffac>0.0) or (xsph_dampfac>0.0 and xsph_stiffac<0.0))
      dserror("Either both or none of the xsph parameters have to be set!");
  }

  // Additional damping force applied in order to achieve static equilibrium solution?
  double damping_factor=DRT::Problem::Instance()->ParticleParams().get<double>("VISCOUS_DAMPING");

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    double density_i = particle_i.density_;

    // determine some useful quantities
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    //******damping (optional)*******************************************************************************************************
    //If required, add artificial viscous damping force in order to determine static equilibrium configurations
    if(damping_factor>0)
    {
      LINALG::Matrix<3,1> dampingforce(true);
      dampingforce.Update(-damping_factor/particle_i.mass_, particle_i.vel_);
      LINALG::Assemble(*accn, dampingforce, particle_i.lm_, particle_i.owner_);
    }

    LINALG::Matrix<3,1> sumj_accn_ij(true);
    LINALG::Matrix<3,1> trvl_sumj_accn_ij(true);

    // loop over the interaction particle list
    boost::unordered_map<int, InterDataPvP>::const_iterator jj;
    for (jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //Boundary particles are not considered here
      if(particle_j.boundarydata_.boundaryparticle_==true)
        continue;

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      double density_j = particle_j.density_;

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);

      const double rho2 = density_i * density_j;
      const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

      LINALG::Matrix<3,1> momentum_ij(true);
      LINALG::Matrix<3,1> p0_momentum_ij(true);
      LINALG::Matrix<3,1> xsph_momentum_ij(true);
      LINALG::Matrix<3,1> momentum_ij_A(true);

      //*******pressure*********************************************************************************************************
      // compute the pressure term
      //(based on Monaghan2005, Eq(3.18); this expression differs from the variant in Espanol2003, Eq(22) in the general case m_i \neq m_j.
      //However, based on a similar derivation as in Espanol2003, Eq(18-22) it has been verified that also this variant can guarantee
      //for energy conservation of the dissipation-free, time-continous problem, if Monaghan2005, Eq(2.18) is applied for density
      // integration [and not Eq(2.17)!])
      const double rhoSquare_j = std::pow(density_j,2);
      const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
      momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);

      //*******transport velocity (optional)***********************************************************************************
      if(background_pressure>0)
      {
        //background pressure term for modified convection velocities/accelerations
        const double p0_gradP0_Rho2 = background_pressure*(1.0/rhoSquare_i + 1.0/rhoSquare_j);  // p_0(1.0/\rho_i^2 + 1.0/\rho_j^2)
        p0_momentum_ij.Update(- p0_gradP0_Rho2, rRelVersor_ij, 1.0);
      }

      if(xsph_dampfac>0 or xsph_stiffac>0)
      {
        //xsph contributions to transport velocity
        double w_ij_damp=interData_ij.w_ij_;
        double w_ij_stiff=weightFunctionHandler_->W(interData_ij.rRelNorm2_, particle_i.radius_/PARTICLE_P0ZHANG);
        xsph_momentum_ij.Update(-xsph_dampfac*w_ij_damp,vRel_ij,1.0);
        xsph_momentum_ij.Update(xsph_stiffac*w_ij_stiff,rRelVersor_ij,1.0);
      }

      if(transport_velocity)
      {
        //additional term div(A) in Navier-Stokes equation due to non-material particle convection
        //with ternsor A_{ab}:=\rho*vel_a*(\tilde{vel}_b-vel_b) according to Adami 2013, Eq. (4).
        LINALG::Matrix<3,3> p0_A_Rho2(true);
        LINALG::Matrix<3,1> p0_A_Rho2_e_ij(true);
        for(int a=0;a<3;a++)
        {
          for(int b=0;b<3;b++)
          {
            p0_A_Rho2(a,b)+=particle_i.vel_(a)*(particle_i.velConv_(b)-particle_i.vel_(b))/density_i;
            p0_A_Rho2(a,b)+=particle_j.vel_(a)*(particle_j.velConv_(b)-particle_j.vel_(b))/density_j;
          }
        }
        p0_A_Rho2_e_ij.MultiplyNN(p0_A_Rho2,rRelVersor_ij);
        if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"NO_VELDIFF_TERM")==false)
        {
          momentum_ij.Update(1.0, p0_A_Rho2_e_ij, 1.0);
        }
        momentum_ij_A.Update(1.0, p0_A_Rho2_e_ij, 1.0);
      }

      //*******viscous forces***************************************************************************************************
      // The following two viscous terms are taken from Espanol2003, Eq(30) and re-expressed in terms of mass densities in an
      // equivalent manner. This is necessary since Espanol2003 only specifies the case m_i=m_j=m.
      // diffusion
      const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
      momentum_ij.Update(dC_rho2rRelNorm2, vRel_ij, 1.0);
      // convection
      const double cCrRelVersorDotVrel_rho2rRelNorm2 =  convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
      momentum_ij.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor_ij, 1.0);

      //*******artificial viscous forces****************************************************************************************
      // The following term represents the forces resulting from artificial viscosity as applied in Adami et al. 2012, Eq. (11).
      // It is only applied for fluid-fluid particle interaction, not for fluid-boundary particle interaction!
      if(artificialViscosity_>0)
      {
        const double h_ij= 0.5*(weightFunctionHandler_->SmoothingLength(particle_i.radius_)+weightFunctionHandler_->SmoothingLength(particle_j.radius_));
        //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
        //TODO: So far, all fluid particles have the same speed of sound. In cases, where each particle has an individual speed of sound,
        // c_ij has to be the average speed of sound of the two particles i and j!
        const double c_ij=particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
        const double density_ij=0.5*(density_i+density_j);
        //The parameter epsilon avoids division by zero when particles get to close. The value epsilon=0.01 is a typical choice (see Adami et al. 2012).
        const double epsilon=0.01;
        const double art_visc_fac=artificialViscosity_*h_ij*c_ij*rRelVersorDotVrel*rRelNorm2 / (density_ij*(rRelNorm2*rRelNorm2+epsilon*h_ij*h_ij));
        momentum_ij.Update(art_visc_fac, rRelVersor_ij, 1.0);
      }

      //*******surface tension**************************************************************************************************
      // The following term applies pairwise interaction forces to model surface tension following Tartakovsky et al. 2016
      if(surfaceTension_>0.0)
      {
        double radius_ij = 0.5*(particle_i.radius_+particle_j.radius_);
        double lambda = 0.0;
        double potential = 0.0;
        SurfTensionInterPot(radius_ij, rRelNorm2, lambda, potential);

        if(potential != 0.0)
        {
          // averaged number density (n_i = \rho_i/m_i)
          double n_ij = 0.5*((density_i/particle_i.mass_)+(density_j/particle_j.mass_));

          double timefac = 1.0;
          if(time >= 0.0)
            timefac = SurfTensionTimeFac(time);

          // pair-wise interaction force
          const double s_ff = 0.5*std::pow(n_ij,-2.0)*(surfaceTension_/lambda);
          LINALG::Matrix<3,1> pf_ij(true);
          pf_ij.Update(-s_ff*potential*timefac, rRelVersor_ij, 1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ij;
          accn_ij.Update(1.0/particle_i.mass_, pf_ij);
          accn_ij.Scale(extMulti);
          sumj_accn_ij.Update(1.0,accn_ij,1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ji;
          accn_ji.Update(1.0/particle_j.mass_, pf_ij);
          accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        accn_ij.Scale(extMulti);
        sumj_accn_ij.Update(1.0,accn_ij,1.0);

        LINALG::Matrix<3,1> accn_ij_A;
        accn_ij_A.Update(generalCoeff_ij, momentum_ij_A);
        accn_ij_A.Scale(extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ij_A, particle_i.lm_, particle_i.owner_);

        if(background_pressure>0)
        {
          LINALG::Matrix<3,1> p0_accn_ij;
          if(freeSurfaceType_==INPAR::PARTICLE::FreeSurface_None)
          {
            p0_accn_ij.Update(generalCoeff_ij, p0_momentum_ij);
            p0_accn_ij.Scale(extMulti);
          }
          else
          {
            double mod_background_pressure=background_pressure;
            if(abs(10.0*particle_i.pressure_)<background_pressure)
              mod_background_pressure=abs(10.0*particle_i.pressure_);

            double dw_ij_tilde=weightFunctionHandler_->DW(interData_ij.rRelNorm2_, particle_i.radius_/PARTICLE_P0ZHANG);
            const double p0_gradP0_Rho2_ij = mod_background_pressure*particle_j.mass_*dw_ij_tilde/rhoSquare_i;
            p0_accn_ij.Update(- p0_gradP0_Rho2_ij, rRelVersor_ij, 1.0);
          }
          trvl_sumj_accn_ij.Update(1.0,p0_accn_ij,1.0);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ij;
          xsph_accn_ij.Update(particle_j.mass_/(0.5*(particle_i.density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          trvl_sumj_accn_ij.Update(1.0,xsph_accn_ij,1.0);
        }
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.dw_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        accn_ji.Scale(-extMulti); // actio = - reactio
        LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);

        LINALG::Matrix<3,1> accn_ji_A;
        accn_ji_A.Update(generalCoeff_ji, momentum_ij_A);
        accn_ji_A.Scale(-extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ji_A, particle_j.lm_, particle_j.owner_);

        if(background_pressure>0)
        {
          LINALG::Matrix<3,1> p0_accn_ji;

          if(freeSurfaceType_==INPAR::PARTICLE::FreeSurface_None)
          {
            p0_accn_ji.Update(generalCoeff_ji, p0_momentum_ij);
            p0_accn_ji.Scale(-extMulti); // actio = - reactio
          }
          else
          {
            double mod_background_pressure=background_pressure;
            if(abs(10.0*particle_j.pressure_)<background_pressure)
              mod_background_pressure=abs(10.0*particle_j.pressure_);

            double dw_ji_tilde=weightFunctionHandler_->DW(interData_ij.rRelNorm2_, particle_j.radius_/PARTICLE_P0ZHANG);
            const double p0_gradP0_Rho2_ji = mod_background_pressure*particle_i.mass_*dw_ji_tilde/rhoSquare_j;
            p0_accn_ji.Update(p0_gradP0_Rho2_ji, rRelVersor_ij, 1.0);
          }
          LINALG::Assemble(*trvl_acc, p0_accn_ji, particle_j.lm_, particle_j.owner_);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ji;
          xsph_accn_ji.Update(-particle_i.mass_/(0.5*(particle_i.density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          LINALG::Assemble(*trvl_acc, xsph_accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }
    }//sum over j
    //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*accn, sumj_accn_ij, particle_i.lm_, particle_i.owner_);
    if(background_pressure>0)
      LINALG::Assemble(*trvl_acc, trvl_sumj_accn_ij, particle_i.lm_, particle_i.owner_);
  }//sum over i

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    Inter_bpvp_acc_var1(accn, trvl_acc, acc_A, time, extMulti);
}


/*----------------------------------------------------------------------*
 | evaluate acceleration - bpvp                            meier 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_bpvp_acc_var1(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> trvl_acc,
    const Teuchos::RCP<Epetra_Vector> acc_A,
    const double time,
    const double extMulti)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }

  // Modified particle transport velocities either due to constant background pressure according to Adami2013 or due to XSPH according to Monaghan1989
  // The terms associated with the constant background pressure are similar to the real pressure terms
  // but with the pressures p_i and p_j replaced by the background pressure. All the related variables
  // and terms are marked by a prefix p0, e.g. p0_momentum_ij etc.. All terms related to XSPH are marked by a prefix xsph
  double background_pressure=-1.0;
  double xsph_dampfac=-1.0;
  double xsph_stiffac=-1.0;
  bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");
  if (trvl_acc != Teuchos::null)
  {
    background_pressure=DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE");
    xsph_dampfac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_DAMPFAC");
    xsph_stiffac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_STIFFAC");

    if(background_pressure<0.0 and xsph_dampfac<0.0 and xsph_stiffac<0.0)
      dserror("Negative values of background pressure, xsph_dampfac and xsph_stifffac! Why is the vector accmod handed in for this case? Set flag TRANSPORT_VELOCITY to no!");

    if(transport_velocity==false)
      dserror("Vector accmod handed in even though TRANSPORT_VELOCITY is set to no?!?");


    if((xsph_dampfac<0.0 and xsph_stiffac>0.0) or (xsph_dampfac>0.0 and xsph_stiffac<0.0))
      dserror("Either both or none of the xsph parameters have to be set!");
  }

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;
    ParticleMF *particle_i = boundaryparticles_[id_i];

    //-1000 means there are no relevant fluid particle neighbors for this boundary particle
    if(particle_i->boundarydata_.pressureMod_==-1000)
      continue;

    //Here, the initial density of the boundary particles is applied
    //(the density of boundary particles remains unchanged since no data is assembled during the simulation to boundary particle DoFs)
    //In Ref. Adami2012, Eq. (28), an alternative density is proposed for calculating a density weighted pressure average.
    //Since we do not apply a density weighted pressure average, this alternative density definition is not required.
    double density_i = particle_i->boundarydata_.densityMod_;

    // determine some useful quantities
    const double rhoSquare_i = std::pow(density_i,2);
    //The pressure terms requires a modified pressure of the boundary particles according to Ref. Adami2012, Eq. (27)
    const double p_Rho2_i = particle_i->boundarydata_.pressureMod_/rhoSquare_i;

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {
      int id_j = jj->first;
      int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      //Atenttion:  particle_i is allways a boundary particle and particle_j is allways a fluid particle.
      //            Forces are only assembled to the fluid particles particle_j!
      //            Everything else is similar to Inter_bpvp_acc(), but with modified boundary particle pressurs and velocities.

      //fluid particle density;
      double density_j = particle_j.density_;

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      //The viscous terms require modified velocities of the boundary particles according to Ref. Adami2012, Eq. (23)
      vRel_ij.Update(1.0, particle_i->boundarydata_.velModVisc_, -1.0, particle_j.vel_);

      LINALG::Matrix<3,1> momentum_ij;
      LINALG::Matrix<3,1> p0_momentum_ij(true);
      LINALG::Matrix<3,1> xsph_momentum_ij(true);
      LINALG::Matrix<3,1> momentum_ij_A;

      //*******pressure*********************************************************************************************************
      //(based on Monaghan2005, Eq(3.8); this expression differs from the variant in Espanol2003, Eq(22) in the general case m_i \neq m_j.
      //However, based on a similar derivation as in Espanol2003, Eq(18-22) it has been verified that also this variant can guarantee
      //for energy conservation of the dissipation-free, time-continous problem, if Monaghan2005, Eq(2.18) is applied for density integration [and not Eq(2.17)!])
      const double rhoSquare_j = std::pow(density_j,2);
      const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
      momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);

      //*******transport velocity (optional)***********************************************************************************
      if(background_pressure>0)
      {
        //background pressure term for modified convection velocities/accelerations
        const double p0_gradP0_Rho2 = background_pressure*(1.0/rhoSquare_i + 1.0/rhoSquare_j);  // p_0(1.0/\rho_i^2 + 1.0/\rho_j^2)
        p0_momentum_ij.Update(- p0_gradP0_Rho2, rRelVersor_ij, 1.0);
      }

      if(xsph_dampfac>0 or xsph_stiffac>0)
      {
        //xsph contributions to transport velocity
        double w_ij_damp=interData_ij.w_ij_;
        double w_ij_stiff=weightFunctionHandler_->W(interData_ij.rRelNorm2_, particle_i->radius_/PARTICLE_P0ZHANG);
        xsph_momentum_ij.Update(-xsph_dampfac*w_ij_damp,vRel_ij,1.0);
        xsph_momentum_ij.Update(xsph_stiffac*w_ij_stiff,rRelVersor_ij,1.0);
      }

      if(transport_velocity)
      {
        //additional term div(A) in Navier-Stokes equation due to non-material particle convection
        //with ternsor A_{ab}:=\rho*vel_a*(\tilde{vel}_b-vel_b) according to Adami 2013, Eq. (4).
        LINALG::Matrix<3,3> p0_A_Rho2(true);
        LINALG::Matrix<3,1> p0_A_Rho2_e_ij(true);
        for(int a=0;a<3;a++)
        {
          for(int b=0;b<3;b++)
          {
            // In the original formulation according to Adami et al.2013, the contributions of boundary particles to the tensor A are simply set to zero
            //(according to an email by Stefan Adami, not mentioned in the paper).
            p0_A_Rho2(a,b)+=particle_j.vel_(a)*(particle_j.velConv_(b)-particle_j.vel_(b))/density_j;
          }
        }
        p0_A_Rho2_e_ij.MultiplyNN(p0_A_Rho2,rRelVersor_ij);
        if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"NO_VELDIFF_TERM")==false)
        {
          momentum_ij.Update(1.0, p0_A_Rho2_e_ij, 1.0);
        }
        momentum_ij_A.Update(1.0, p0_A_Rho2_e_ij, 1.0);
      }

      if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip)
      {
        //*******viscous forces***************************************************************************************************
        //Viscous interaction with boundary particles is only required in case of a no-slip boundary condition!
        const double rRelNorm2 = interData_ij.rRelNorm2_;
        const double rho2 = density_i * density_j;
        const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);
        // The following two viscous terms are taken from Espanol2003, Eq(30) and re-expressed in terms of mass densities in an equivalent manner.
        // This is necessary since Espanol2003 only specifies the case m_i=m_j=m.
        // diffusion
        const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
        momentum_ij.Update(dC_rho2rRelNorm2, vRel_ij, 1.0);
        // convection
        const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
        momentum_ij.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor_ij, 1.0);

        //*******artificial viscous forces****************************************************************************************
        // The following term represents the forces resulting from artificial viscosity as applied in Adami et al. 2012, Eq. (11).
        // It is only applied for fluid-fluid particle interaction, not for fluid-boundary particle interaction!
        if(artificialViscosity_>0)
        {
          const double h_ij= 0.5*(weightFunctionHandler_->SmoothingLength(particle_i->radius_)+weightFunctionHandler_->SmoothingLength(particle_j.radius_));
          //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
          //TODO: So far, all fluid particles have the same speed of sound. In cases, where each particle has an individual speed of sound,
          // c_ij has to be the average speed of sound of the two particles i and j!
          const double c_ij=particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
          const double density_ij=0.5*(density_i+density_j);
          //The parameter epsilon avoids division by zero when particles get to close. The value epsilon=0.01 is a typical choice (see Adami et al. 2012).
          const double epsilon=0.01;
          const double art_visc_fac=artificialViscosity_*h_ij*c_ij*rRelVersorDotVrel*rRelNorm2 / (density_ij*(rRelNorm2*rRelNorm2+epsilon*h_ij*h_ij));
          momentum_ij.Update(art_visc_fac, rRelVersor_ij, 1.0);
        }
      }

      //*******surface tension**************************************************************************************************
      // The following term applies pairwise interaction forces to model surface tension following Tartakovsky et al. 2016
      if(surfaceTension_>0.0)
      {
        double radius_ij = 0.5*(particle_i->radius_+particle_j.radius_);
        double lambda = 0.0;
        double potential = 0.0;
        SurfTensionInterPot(radius_ij, interData_ij.rRelNorm2_, lambda, potential);

        if(potential != 0.0)
        {
          // averaged number density (n_i = \rho_i/m_i)
          double n_ij = 0.5*((density_i/particle_i->mass_)+(density_j/particle_j.mass_));

          // convert static contact anglee in radians
          double theta_0 = staticContactAngle_*PI/180.0;

          double timefac = 1.0;
          if(time >= 0.0)
            timefac = SurfTensionTimeFac(time);

          // pair-wise interaction force
          const double s_sf = 0.5*std::pow(n_ij,-2.0)*(surfaceTension_/lambda)*(1+0.5*std::cos(theta_0));
          LINALG::Matrix<3,1> pf_ij(true);
          pf_ij.Update(-s_sf*potential*timefac, rRelVersor_ij, 1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ji;
          accn_ji.Update(1.0/particle_j.mass_, pf_ij);
          accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.dw_ji_ * particle_i->mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        accn_ji.Scale(-extMulti); // actio = - reactio
        LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);

        LINALG::Matrix<3,1> accn_ji_A;
        accn_ji_A.Update(generalCoeff_ji, momentum_ij_A);
        accn_ji_A.Scale(-extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ji_A, particle_j.lm_, particle_j.owner_);

        if(background_pressure>0)
        {
          LINALG::Matrix<3,1> p0_accn_ji;
          if(freeSurfaceType_==INPAR::PARTICLE::FreeSurface_None)
          {
            p0_accn_ji.Update(generalCoeff_ji, p0_momentum_ij);
            p0_accn_ji.Scale(-extMulti); // actio = - reactio
          }
          else
          {
            double mod_background_pressure=background_pressure;
            if(abs(10.0*particle_j.pressure_)<background_pressure)
              mod_background_pressure=abs(10.0*particle_j.pressure_);

            double dw_ji_tilde=weightFunctionHandler_->DW(interData_ij.rRelNorm2_, particle_j.radius_/PARTICLE_P0ZHANG);
            const double p0_gradP0_Rho2_ji = mod_background_pressure*particle_i->mass_*dw_ji_tilde/rhoSquare_j;
            p0_accn_ji.Update(p0_gradP0_Rho2_ji, rRelVersor_ij, 1.0);
          }
          LINALG::Assemble(*trvl_acc, p0_accn_ji, particle_j.lm_, particle_j.owner_);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ji;
          xsph_accn_ji.Update(-particle_i->mass_/(0.5*(particle_i->density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          LINALG::Assemble(*trvl_acc, xsph_accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }
    }//loop over j
  }//loop over i
}


/*------------------------------------------------------------------------------------------*
 | evaluate acceleration - pvp: alternative variant according to Adami2013     meier 05/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_acc_var2(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> trvl_acc,
    const Teuchos::RCP<Epetra_Vector> acc_A,
    const double time,
    const double extMulti)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_acc");

  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }

  // Modified particle transport velocities either due to constant background pressure according to Adami2013 or due to XSPH according to Monaghan1989
  // The terms associated with the constant background pressure are similar to the real pressure terms
  // but with the pressures p_i and p_j replaced by the background pressure. All the related variables
  // and terms are marked by a prefix p0, e.g. p0_momentum_ij etc.. All terms related to XSPH are marked by a prefix xsph
  double background_pressure=-1.0;
  double xsph_dampfac=-1.0;
  double xsph_stiffac=-1.0;
  bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");
  if (trvl_acc != Teuchos::null)
  {
    background_pressure=DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE");
    xsph_dampfac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_DAMPFAC");
    xsph_stiffac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_STIFFAC");

    if(background_pressure<0.0 and xsph_dampfac<0.0 and xsph_stiffac<0.0)
      dserror("Negative values of background pressure, xsph_dampfac and xsph_stifffac! Why is the vector accmod handed in for this case? Set flag TRANSPORT_VELOCITY to no!");

    if(transport_velocity==false)
      dserror("Vector accmod handed in even though TRANSPORT_VELOCITY is set to no?!?");

    if((xsph_dampfac<0.0 and xsph_stiffac>0.0) or (xsph_dampfac>0.0 and xsph_stiffac<0.0))
      dserror("Either both or none of the xsph parameters have to be set!");
  }

  // Additional damping force applied in order to achieve static equilibrium solution?
  const double damping_factor=DRT::Problem::Instance()->ParticleParams().get<double>("VISCOUS_DAMPING");
  const double viscosity=particle_algorithm_->ExtParticleMat()->dynamicViscosity_;

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    double density_i = particle_i.density_;
    double pressure_i = particle_i.pressure_;
    double V_i=particle_i.mass_/density_i;

    //******damping (optional)*******************************************************************************************************
    //If required, add artificial viscous damping force in order to determine static equilibrium configurations
    if(damping_factor>0)
    {
      LINALG::Matrix<3,1> dampingforce(true);
      dampingforce.Update(-damping_factor/particle_i.mass_, particle_i.vel_);
      LINALG::Assemble(*accn, dampingforce, particle_i.lm_, particle_i.owner_);
    }

    LINALG::Matrix<3,1> sumj_accn_ij(true);
    LINALG::Matrix<3,1> trvl_sumj_accn_ij(true);

    // loop over the interaction particle list
    boost::unordered_map<int, InterDataPvP>::const_iterator jj;
    for (jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //Boundary particles are not considered here
      if(particle_j.boundarydata_.boundaryparticle_==true)
        continue;

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      double density_j = particle_j.density_;
      double pressure_j = particle_j.pressure_;
      double V_j=particle_j.mass_/density_j;

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);

      LINALG::Matrix<3,1> momentum_ij(true);
      LINALG::Matrix<3,1> p0_momentum_ij(true);
      LINALG::Matrix<3,1> xsph_momentum_ij(true);
      LINALG::Matrix<3,1> momentum_ij_A(true);
      LINALG::Matrix<3,1> momentum_artvisc_ij(true);

      //*******pressure*********************************************************************************************************
      // compute the pressure term

      const double tilde_p_ij = (density_i*pressure_j+density_j*pressure_i)/(density_i+density_j);
      momentum_ij.Update(- tilde_p_ij, rRelVersor_ij, 1.0);

      //*******transport velocity (optional)***********************************************************************************
      if(background_pressure>0)
      {
        //background pressure term for modified convection velocities/accelerations
        p0_momentum_ij.Update(- background_pressure, rRelVersor_ij, 1.0);
      }

      if(xsph_dampfac>0 or xsph_stiffac>0)
      {
        //xsph contributions to transport velocity
        double w_ij_damp=interData_ij.w_ij_;
        double w_ij_stiff=weightFunctionHandler_->W(interData_ij.rRelNorm2_, particle_i.radius_/PARTICLE_P0ZHANG);
        xsph_momentum_ij.Update(-xsph_dampfac*w_ij_damp,vRel_ij,1.0);
        xsph_momentum_ij.Update(xsph_stiffac*w_ij_stiff,rRelVersor_ij,1.0);
      }

      if(transport_velocity)
      {
        //additional term div(A) in Navier-Stokes equation due to non-material particle convection
        //with ternsor A_{ab}:=\rho*vel_a*(\tilde{vel}_b-vel_b) according to Adami 2013, Eq. (4).
        LINALG::Matrix<3,3> p0_A(true);
        LINALG::Matrix<3,1> p0_A_e_ij(true);
        for(int a=0;a<3;a++)
        {
          for(int b=0;b<3;b++)
          {
            p0_A(a,b)+=density_i*particle_i.vel_(a)*(particle_i.velConv_(b)-particle_i.vel_(b));
            p0_A(a,b)+=density_j*particle_j.vel_(a)*(particle_j.velConv_(b)-particle_j.vel_(b));
          }
        }
        p0_A_e_ij.MultiplyNN(p0_A,rRelVersor_ij);

        if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"NO_VELDIFF_TERM")==false)
        {
          momentum_ij.Update(0.5, p0_A_e_ij, 1.0);
        }
        momentum_ij_A.Update(0.5, p0_A_e_ij, 1.0);

      }

      //*******viscous forces***************************************************************************************************
      momentum_ij.Update(viscosity/rRelNorm2, vRel_ij, 1.0);

      //*******artificial viscous forces****************************************************************************************
      // The following term represents the forces resulting from artificial viscosity as applied in Adami et al. 2012, Eq. (11).
      // It is only applied for fluid-fluid particle interaction, not for fluid-boundary particle interaction!
      if(artificialViscosity_>0)
      {
        const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

        const double h_ij= 0.5*(weightFunctionHandler_->SmoothingLength(particle_i.radius_)+weightFunctionHandler_->SmoothingLength(particle_j.radius_));
        //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
        //TODO: So far, all fluid particles have the same speed of sound. In cases, where each particle has an individual speed of sound,
        // c_ij has to be the average speed of sound of the two particles i and j!
        const double c_ij=particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
        const double density_ij=0.5*(density_i+density_j);
        //The parameter epsilon avoids division by zero when particles get to close. The value epsilon=0.01 is a typical choice (see Adami et al. 2012).
        const double epsilon=0.01;
        //Attention: Here, we use an extra data container for the momentum contribution from artificial viscosity (which was not required in Inter_pvp_acc_var1(...))
        //since this term needs a different pre-factor than the other momentum contributions.
        const double art_visc_fac=artificialViscosity_*h_ij*c_ij*rRelVersorDotVrel*rRelNorm2 / (density_ij*(rRelNorm2*rRelNorm2+epsilon*h_ij*h_ij));
        momentum_artvisc_ij.Update(art_visc_fac, rRelVersor_ij, 1.0);
      }

      //*******surface tension**************************************************************************************************
      // The following term applies pairwise interaction forces to model surface tension following Tartakovsky et al. 2016
      if(surfaceTension_>0.0)
      {
        double radius_ij = 0.5*(particle_i.radius_+particle_j.radius_);
        double lambda = 0.0;
        double potential = 0.0;
        SurfTensionInterPot(radius_ij, rRelNorm2, lambda, potential);

        if(potential != 0.0)
        {
          // averaged number density (n_i = \rho_i/m_i)
          double n_ij = 0.5*((density_i/particle_i.mass_)+(density_j/particle_j.mass_));

          double timefac = 1.0;
          if(time >= 0.0)
            timefac = SurfTensionTimeFac(time);

          // pair-wise interaction force
          const double s_ff = 0.5*std::pow(n_ij,-2.0)*(surfaceTension_/lambda);
          LINALG::Matrix<3,1> pf_ij(true);
          pf_ij.Update(-s_ff*potential*timefac, rRelVersor_ij, 1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ij;
          accn_ij.Update(1.0/particle_i.mass_, pf_ij);
          accn_ij.Scale(extMulti);
          sumj_accn_ij.Update(1.0,accn_ij,1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ji;
          accn_ji.Update(1.0/particle_j.mass_, pf_ij);
          accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = (V_i*V_i+V_j*V_j)*interData_ij.dw_ij_ / particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        if(artificialViscosity_>0)
        {
          const double generalCoeff_artvisc_ij = particle_j.mass_*interData_ij.dw_ij_;
          accn_ij.Update(generalCoeff_artvisc_ij, momentum_artvisc_ij, 1.0);
        }
        accn_ij.Scale(extMulti);
        sumj_accn_ij.Update(1.0,accn_ij,1.0);

        LINALG::Matrix<3,1> accn_ij_A;
        accn_ij_A.Update(generalCoeff_ij, momentum_ij_A);
        accn_ij_A.Scale(extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ij_A, particle_i.lm_, particle_i.owner_);

        if(background_pressure>0)
        {
          if(freeSurfaceType_!=INPAR::PARTICLE::FreeSurface_None)
            dserror("Background pressure not available for free-surface flow in var2 so far!");

          LINALG::Matrix<3,1> p0_accn_ij;
          p0_accn_ij.Update(generalCoeff_ij, p0_momentum_ij);
          p0_accn_ij.Scale(extMulti);
          trvl_sumj_accn_ij.Update(1.0,p0_accn_ij,1.0);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ij;
          xsph_accn_ij.Update(particle_j.mass_/(0.5*(particle_i.density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          trvl_sumj_accn_ij.Update(1.0,xsph_accn_ij,1.0);
        }
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = (V_i*V_i+V_j*V_j)*interData_ij.dw_ji_ / particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        if(artificialViscosity_>0)
        {
          const double generalCoeff_artvisc_ji = particle_i.mass_*interData_ij.dw_ji_;
          accn_ji.Update(generalCoeff_artvisc_ji, momentum_artvisc_ij, 1.0);
        }
        accn_ji.Scale(-extMulti); // actio = - reactio
        LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);

        LINALG::Matrix<3,1> accn_ji_A;
        accn_ji_A.Update(generalCoeff_ji, momentum_ij_A);
        accn_ji_A.Scale(-extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ji_A, particle_j.lm_, particle_j.owner_);

        if(background_pressure>0)
        {
          if(freeSurfaceType_!=INPAR::PARTICLE::FreeSurface_None)
            dserror("Background pressure not available for free-surface flow in var2 so far!");

          LINALG::Matrix<3,1> p0_accn_ji;
          p0_accn_ji.Update(generalCoeff_ji, p0_momentum_ij);
          p0_accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*trvl_acc, p0_accn_ji, particle_j.lm_, particle_j.owner_);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ji;
          xsph_accn_ji.Update(-particle_i.mass_/(0.5*(particle_i.density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          LINALG::Assemble(*trvl_acc, xsph_accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }
    }//sum over j
    //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*accn, sumj_accn_ij, particle_i.lm_, particle_i.owner_);
    if(background_pressure>0)
      LINALG::Assemble(*trvl_acc, trvl_sumj_accn_ij, particle_i.lm_, particle_i.owner_);
  }//sum over i

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    Inter_bpvp_acc_var2(accn, trvl_acc, acc_A, extMulti);
}

/*-----------------------------------------------------------------------------------------*
 | evaluate acceleration - bpvp: alternative variant according to Adami2013   meier 02/17  |
 *-----------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_bpvp_acc_var2(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> trvl_acc,
    const Teuchos::RCP<Epetra_Vector> acc_A,
    const double time,
    const double extMulti)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }

  // Modified particle transport velocities either due to constant background pressure according to Adami2013 or due to XSPH according to Monaghan1989
  // The terms associated with the constant background pressure are similar to the real pressure terms
  // but with the pressures p_i and p_j replaced by the background pressure. All the related variables
  // and terms are marked by a prefix p0, e.g. p0_momentum_ij etc.. All terms related to XSPH are marked by a prefix xsph
  double background_pressure=-1.0;
  double xsph_dampfac=-1.0;
  double xsph_stiffac=-1.0;
  bool transport_velocity=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY");
  if (trvl_acc != Teuchos::null)
  {
    background_pressure=DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE");
    xsph_dampfac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_DAMPFAC");
    xsph_stiffac=DRT::Problem::Instance()->ParticleParams().get<double>("XSPH_STIFFAC");

    if(background_pressure<0.0 and xsph_dampfac<0.0 and xsph_stiffac<0.0)
      dserror("Negative values of background pressure, xsph_dampfac and xsph_stifffac! Why is the vector accmod handed in for this case? Set flag TRANSPORT_VELOCITY to no!");

    if(transport_velocity==false)
      dserror("Vector accmod handed in even though TRANSPORT_VELOCITY is set to no?!?");

    if((xsph_dampfac<0.0 and xsph_stiffac>0.0) or (xsph_dampfac>0.0 and xsph_stiffac<0.0))
      dserror("Either both or none of the xsph parameters have to be set!");
  }

  const double viscosity=particle_algorithm_->ExtParticleMat()->dynamicViscosity_;

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;
    ParticleMF *particle_i = boundaryparticles_[id_i];

    //-1000 means there are no relevant fluid particle neighbors for this boundary particle
    if(particle_i->boundarydata_.pressureMod_==-1000)
      continue;

    //Here, the initial density of the boundary particles is applied
    //(the density of boundary particles remains unchanged since no data is assembled during the simulation to boundary particle DoFs)
    //In Ref. Adami2012, Eq. (28), an alternative density is proposed for calculating a density weighted pressure average.
    //Since we do not apply a density weighted pressure average, this alternative density definition is not required.
    //TODO:Check the difference with modified density!!!
    double density_i = particle_i->boundarydata_.densityMod_;
    double pressure_i =  particle_i->boundarydata_.pressureMod_;
    double V_i=particle_i->mass_/density_i;

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {
      int id_j = jj->first;
      int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      //Atenttion:  particle_i is allways a boundary particle and particle_j is allways a fluid particle.
      //            Forces are only assembled to the fluid particles particle_j!
      //            Everything else is similar to Inter_bpvp_acc(), but with modified boundary particle pressurs and velocities.

      double density_j = particle_j.density_;
      double pressure_j = particle_j.pressure_;
      double V_j=particle_j.mass_/density_j;

      // --- extract general data --- //
      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      //The viscous terms require modified velocities of the boundary particles according to Ref. Adami2012, Eq. (23)
      vRel_ij.Update(1.0, particle_i->boundarydata_.velModVisc_, -1.0, particle_j.vel_);

      LINALG::Matrix<3,1> momentum_ij(true);
      LINALG::Matrix<3,1> p0_momentum_ij(true);
      LINALG::Matrix<3,1> xsph_momentum_ij(true);
      LINALG::Matrix<3,1> momentum_ij_A(true);
      LINALG::Matrix<3,1> momentum_artvisc_ij(true);

      //*******pressure*********************************************************************************************************
      // compute the pressure term
      const double tilde_p_ij = (density_i*pressure_j+density_j*pressure_i)/(density_i+density_j);
      momentum_ij.Update(- tilde_p_ij, rRelVersor_ij, 1.0);

      //*******transport velocity (optional)***********************************************************************************
      if(background_pressure>0)
      {
        //background pressure term for modified convection velocities/accelerations
        p0_momentum_ij.Update(- background_pressure, rRelVersor_ij, 1.0);
      }

      if(xsph_dampfac>0 or xsph_stiffac>0)
      {
        //xsph contributions to transport velocity
        double w_ij_damp=interData_ij.w_ij_;
        double w_ij_stiff=weightFunctionHandler_->W(interData_ij.rRelNorm2_, particle_i->radius_/PARTICLE_P0ZHANG);
        xsph_momentum_ij.Update(-xsph_dampfac*w_ij_damp,vRel_ij,1.0);
        xsph_momentum_ij.Update(xsph_stiffac*w_ij_stiff,rRelVersor_ij,1.0);
      }

      if(transport_velocity)
      {
        //additional term div(A) in Navier-Stokes equation due to non-material particle convection
        //with ternsor A_{ab}:=\rho*vel_a*(\tilde{vel}_b-vel_b) according to Adami 2013, Eq. (4).
        LINALG::Matrix<3,3> p0_A(true);
        LINALG::Matrix<3,1> p0_A_e_ij(true);
        for(int a=0;a<3;a++)
        {
          for(int b=0;b<3;b++)
          {
            // In the original formulation according to Adami et al.2013, the contributions of boundary particles to the tensor A are simply set to zero
            //(according to an email by Stefan Adami, not mentioned in the paper).
            p0_A(a,b)+=density_j*particle_j.vel_(a)*(particle_j.velConv_(b)-particle_j.vel_(b));
          }
        }
        p0_A_e_ij.MultiplyNN(p0_A,rRelVersor_ij);

        if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"NO_VELDIFF_TERM")==false)
        {
          momentum_ij.Update(0.5, p0_A_e_ij, 1.0);
        }
        momentum_ij_A.Update(0.5, p0_A_e_ij, 1.0);
      }

      if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip)
      {
        //*******viscous forces***************************************************************************************************
        //Viscous interaction with boundary particles is only required in case of a no-slip boundary condition!

          momentum_ij.Update(viscosity/rRelNorm2, vRel_ij, 1.0);

        //*******artificial viscous forces****************************************************************************************
        // The following term represents the forces resulting from artificial viscosity as applied in Adami et al. 2012, Eq. (11).
        // It is only applied for fluid-fluid particle interaction, not for fluid-boundary particle interaction!
        if(artificialViscosity_>0)
        {
          const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

          const double h_ij= 0.5*(weightFunctionHandler_->SmoothingLength(particle_i->radius_)+weightFunctionHandler_->SmoothingLength(particle_j.radius_));
          //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
          //TODO: So far, all fluid particles have the same speed of sound. In cases, where each particle has an individual speed of sound,
          // c_ij has to be the average speed of sound of the two particles i and j!
          const double c_ij=particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
          const double density_ij=0.5*(density_i+density_j);
          //The parameter epsilon avoids division by zero when particles get to close. The value epsilon=0.01 is a typical choice (see Adami et al. 2012).
          const double epsilon=0.01;
          //Attention: Here, we use an extra data container for the momentum contribution from artificial viscosity (which was not required in Inter_pvp_acc_var1(...))
          //since this term needs a different pre-factor than the other momentum contributions.
          const double art_visc_fac=artificialViscosity_*h_ij*c_ij*rRelVersorDotVrel*rRelNorm2 / (density_ij*(rRelNorm2*rRelNorm2+epsilon*h_ij*h_ij));
          momentum_artvisc_ij.Update(art_visc_fac, rRelVersor_ij, 1.0);
        }
      }

      //*******surface tension**************************************************************************************************
      // The following term applies pairwise interaction forces to model surface tension following Tartakovsky et al. 2016
      if(surfaceTension_>0.0)
      {
        double radius_ij = 0.5*(particle_i->radius_+particle_j.radius_);
        double lambda = 0.0;
        double potential = 0.0;
        SurfTensionInterPot(radius_ij, interData_ij.rRelNorm2_, lambda, potential);

        if(potential != 0.0)
        {
          // averaged number density (n_i = \rho_i/m_i)
          double n_ij = 0.5*((density_i/particle_i->mass_)+(density_j/particle_j.mass_));

          // convert static contact anglee in radians
          double theta_0 = staticContactAngle_*PI/180.0;

          double timefac = 1.0;
          if(time >= 0.0)
            timefac = SurfTensionTimeFac(time);

          // pair-wise interaction force
          const double s_sf = 0.5*std::pow(n_ij,-2.0)*(surfaceTension_/lambda)*(1+0.5*std::cos(theta_0));
          LINALG::Matrix<3,1> pf_ij(true);
          pf_ij.Update(-s_sf*potential*timefac, rRelVersor_ij, 1.0);

          // assemble and write
          LINALG::Matrix<3,1> accn_ji;
          accn_ji.Update(1.0/particle_j.mass_, pf_ij);
          accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        double generalCoeff_ji = 1.0;
        generalCoeff_ji = (V_i*V_i+V_j*V_j)*interData_ij.dw_ji_ / particle_j.mass_;

        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        if(artificialViscosity_>0)
        {
          const double generalCoeff_artvisc_ji = particle_i->mass_*interData_ij.dw_ji_;
          accn_ji.Update(generalCoeff_artvisc_ji, momentum_artvisc_ij, 1.0);
        }
        accn_ji.Scale(-extMulti); // actio = - reactio
        LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);

        LINALG::Matrix<3,1> accn_ji_A;
        accn_ji_A.Update(generalCoeff_ji, momentum_ij_A);
        accn_ji_A.Scale(-extMulti);
        if(acc_A!=Teuchos::null)
          LINALG::Assemble(*acc_A, accn_ji_A, particle_j.lm_, particle_j.owner_);

        if(background_pressure>0)
        {
          if(freeSurfaceType_!=INPAR::PARTICLE::FreeSurface_None)
            dserror("Background pressure not available for free-surface flow in var2 so far!");

          LINALG::Matrix<3,1> p0_accn_ji;
          p0_accn_ji.Update(generalCoeff_ji, p0_momentum_ij);
          p0_accn_ji.Scale(-extMulti); // actio = - reactio
          LINALG::Assemble(*trvl_acc, p0_accn_ji, particle_j.lm_, particle_j.owner_);
        }

        if(xsph_dampfac>0 or xsph_stiffac>0)
        {
          LINALG::Matrix<3,1> xsph_accn_ji;
          xsph_accn_ji.Update(-particle_i->mass_/(0.5*(particle_i->density_+particle_j.density_)), xsph_momentum_ij, 1.0);
          LINALG::Assemble(*trvl_acc, xsph_accn_ji, particle_j.lm_, particle_j.owner_);
        }
      }
    }//loop over j
  }//loop over i
}

/*----------------------------------------------------------------------*
 | Computes Laplace operator for vector v - pvp            meier 03/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_Laplace_x(
    const ParticleMF& particle_j,
    const InterDataPvP& interData_ij,
    const LINALG::Matrix<3,1>& x_ij,
    LINALG::Matrix<3,1>& Laplace_x_ij)
{
  //This method determines the contribution of the particle pair i and j to the Laplace operator applied to an arbitrary
  //vector x according to Espanol2003 Eq. (28, first line). This formula is valid independent from the spatial dimension (1D, 2D or 3D).

  Laplace_x_ij=x_ij;
  double F_ij=interData_ij.dw_ij_/interData_ij.rRelNorm2_;
  double d_j=particle_j.density_/particle_j.mass_;
  Laplace_x_ij.Scale(-2.0*F_ij/d_j);
}

/*-----------------------------------------------------------------------------*
 | Computes Gradient of divergence operator for vector v - pvp    meier 03/17  |
 *-----------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_GradDiv_x(
    const ParticleMF& particle_j,
    const InterDataPvP& interData_ij,
    const LINALG::Matrix<3,1>& x_ij,
    LINALG::Matrix<3,1>& GradDiv_x_ij)
{
  //This method determines the contribution of the particle pair i and j to the Laplace operator applied to an arbitrary
  //vector x according to Espanol2003 Eq. (28, second line). In this reference, the formula is only stated for the 3D
  //case. By conducting the derivation in Appendix A for the 2D as well as for the 1D case, it can be shown that in these
  //cases the original coefficient 5 (3D) has to be set to 4 (2D) or 3 (1D), respectively.

  double dim_fac=0.0;
  switch (WF_DIM_)
  {
    case INPAR::PARTICLE::WF_3D :
    {
      dim_fac=5.0;
      break;
    }
    case INPAR::PARTICLE::WF_2D :
    {
      dim_fac=4.0;
      break;
    }
    case INPAR::PARTICLE::WF_1D :
    {
      dim_fac=3.0;
      break;
    }
    default :
    {
      dserror("Only the problem dimensions 1, 2 and 3 are possible!");
    }
  }

  double eij_dot_xij=x_ij.Dot(interData_ij.rRelVersor_ij_);
  double F_ij=interData_ij.dw_ij_/interData_ij.rRelNorm2_;
  double d_j=particle_j.density_/particle_j.mass_;
  GradDiv_x_ij.Update(-1.0,x_ij,0.0);
  GradDiv_x_ij.Update(dim_fac*eij_dot_xij,interData_ij.rRelVersor_ij_,1.0);
  GradDiv_x_ij.Scale(-F_ij/d_j);
}

/*----------------------------------------------------------------------*
 | evaluate specEnthalpyDot - pvp                          katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_specEnthalpyDot(
    const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
    const double extMulti)
{
  //checks
  if (specEnthalpyDotn == Teuchos::null)
  {
    dserror("specEnthalpyDotn is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      dserror("How did you get here? The combination of boundary particles and thermal problems is not possible yet!");

    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //Boundary particles are not considered here
      if(particle_j.boundarydata_.boundaryparticle_==true)
        dserror("How did you get here? The combination of boundary particles and thermal problems is not possible yet!");

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      const double rho2 = density_i * density_j;
      // assemble and write
      const double deltaT_ij = particle_i.temperature_ - particle_j.temperature_;
      const double divT_rho_ij = 2 * particle_algorithm_->ExtParticleMat()->thermalConductivity_ *
          deltaT_ij / (rho2 * rRelNorm2);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;

        // --- specific enthalpy --- //

        double specEnthalpyDotn_ij = generalCoeff_ij * divT_rho_ij;
        specEnthalpyDotn_ij *= extMulti;
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.dw_ji_ * particle_i.mass_;
        // assemble and write
        double specEnthalpyDotn_ji = generalCoeff_ji * divT_rho_ij;
        specEnthalpyDotn_ji *= -extMulti; // actio = - reactio
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ji, particle_j.gid_, particle_j.owner_);
      }
    }
  }
}

/*--------------------------------------------------------------------------*
 | evaluate specEnthalpyDot - pvhs                             katta 10/16  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvhs_specEnthalpyDot(
    const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
    const double extMulti)
{
//checks
if (specEnthalpyDotn == Teuchos::null)
{
  dserror("specEnthalpyDotn is empty");
}

// loop over the particles (no superpositions)
for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_hs_.size(); ++lidNodeRow_i)
{
  // determine the particle_i
  const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
  const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
  const ParticleMF& particle_i = colParticles_[lidNodeCol_i];
  //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
  const double density_i = particle_i.density_;

  double specEnthalpyDot_i = 0.0;
  // loop over the interaction particle list
  for (boost::unordered_map<int, Teuchos::RCP<HeatSource> >::const_iterator jj = neighbours_hs_[lidNodeRow_i].begin(); jj != neighbours_hs_[lidNodeRow_i].end(); ++jj)
  {
    const Teuchos::RCP<HeatSource> hs = jj->second;
    specEnthalpyDot_i += (hs->QDot_)/density_i;
  }
  specEnthalpyDot_i *= extMulti;
  LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDot_i, particle_i.gid_, particle_i.owner_);
}
}

/*--------------------------------------------------------------------------*
 | compute \sum m * W (usually the density) - mesh free style  katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW(
    const Teuchos::RCP<Epetra_Vector> mW,
    const Teuchos::RCP<Epetra_Vector> unity_vec,
    const double extMulti)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW");

  //checks
  if (mW == Teuchos::null)
  {
    dserror("mW is empty");
  }

  // erase the vectors
  mW->PutScalar(0.0);

  if(unity_vec!=Teuchos::null)
    unity_vec->PutScalar(0.0);

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    // auto-interaction
    double mW_ii = weightFunctionHandler_->W0(particle_i.radius_) * particle_i.mass_;
    mW_ii *= extMulti;
    LINALG::Assemble(*mW, mW_ii, particle_i.gid_, particle_i.owner_);

    double unity_vec_aux=0.0;
    if(unity_vec!=Teuchos::null)
    {
      unity_vec_aux=mW_ii/particle_i.density_;
      LINALG::Assemble(*unity_vec, unity_vec_aux, particle_i.gid_, particle_i.owner_);
    }

    double sumj_mW_ij=0.0;
    double sumj_m_rho_W_ij=0.0;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //Boundary particles are not considered here
      if(particle_j.boundarydata_.boundaryparticle_==true)
        continue;

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.w_ij_ != 0)
      {
        // assemble and write
        double mW_ij = interData_ij.w_ij_ * particle_j.mass_;
        mW_ij *= extMulti;
        sumj_mW_ij+=mW_ij;
        if(unity_vec!=Teuchos::null)
        {
          sumj_m_rho_W_ij+=mW_ij/particle_j.density_;
        }
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.w_ji_ != 0)
      {
        // assemble and write
        double mW_ji = interData_ij.w_ji_ * particle_i.mass_;
        mW_ji *= extMulti;
        LINALG::Assemble(*mW, mW_ji, particle_j.gid_, particle_j.owner_);

        if(unity_vec!=Teuchos::null)
        {
          unity_vec_aux=mW_ji/particle_i.density_;
          LINALG::Assemble(*unity_vec, unity_vec_aux, particle_j.gid_, particle_j.owner_);
        }
      }
    }//loop over j
    //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*mW, sumj_mW_ij, particle_i.gid_, particle_i.owner_);
    if(unity_vec!=Teuchos::null)
      LINALG::Assemble(*unity_vec, sumj_m_rho_W_ij, particle_i.gid_, particle_i.owner_);
  }//loop over i

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    MF_mW_Boundary(mW,unity_vec,extMulti);

  return;
}

/*------------------------------------------------------------------------------------------*
 | compute \sum m * W (usually the density) - boundary particle contributions  meier 02/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW_Boundary(
    const Teuchos::RCP<Epetra_Vector> mW,
    const Teuchos::RCP<Epetra_Vector> unity_vec,
    const double extMulti)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW_Boundary");

  //checks
  if (mW == Teuchos::null)
  {
    dserror("mW is empty");
  }

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;
    ParticleMF *particle_i = boundaryparticles_[id_i];

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {
      int id_j = jj->first;
      int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      //Atenttion:  particle_i is allways a boundary particle and particle_j is allways a fluid particle.
      //            Forces are only assembled to the fluid particles particle_j!
      //            Everything else is very similar to MF_mW().

      // write on particle j if appropriate specializing the quantities
      if ( interData_ij.dw_ji_ != 0)
      {
        // assemble and write
        double mW_ji = interData_ij.w_ji_ * particle_i->mass_;
        mW_ji *= extMulti;
        LINALG::Assemble(*mW, mW_ji, particle_j.gid_, particle_j.owner_);
        if(unity_vec!=Teuchos::null)
        {
          double unity_vec_aux=0.0;
#ifndef PARTICLE_BOUNDARYDENSITY
          unity_vec_aux=mW_ji/initDensity_;
#else
          unity_vec_aux=mW_ji/particle_i->boundarydata_.densityMod_;
#endif


          LINALG::Assemble(*unity_vec, unity_vec_aux, particle_j.gid_, particle_j.owner_);
        }
      }
    }
  }

  // loop over all boundary particles (not only the one with fluid neighbors!) and set density and color field
  // loop over the boundary particles
  for (std::map<int,ParticleMF* >::iterator ii = boundaryparticles_.begin(); ii!=boundaryparticles_.end(); ++ii)
  {
    ParticleMF* particle_i = ii->second;

    //Since the density of the boundary particles particle_i should remain unchanged, we have to assemble the initial density
    //value again after all DoFs of the vector mW have been cleared above. This procedure is not required for boundary particle density calculation
    //(or position update) via time integration since vanishing entries for densityDot in the boundary particle DoFs automatically yield an unchanged density.
    double density = initDensity_;
    LINALG::Assemble(*mW, density, particle_i->gid_, particle_i->owner_);
    if(unity_vec!=Teuchos::null)
    {
      double unity = 1.0;
      LINALG::Assemble(*unity_vec, unity, particle_i->gid_, particle_i->owner_);
    }
  }

  return;
}

/*------------------------------------------------------------------------------------*
 | re-initialize density in case of free-surface flow - mesh free style  meier 09/17  |
 *------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_ReInitDensity(
    const Teuchos::RCP<Epetra_Vector> density,
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::MF_ReInitDensity");

  //checks
  if (density == Teuchos::null)
  {
    dserror("density empty");
  }

  // erase the vectors
  density->PutScalar(0.0);

  if(freeSurfaceType==INPAR::PARTICLE::InteriorReinitialization)
  {
    //loop over particles and re-initialize interior particles (FS_NONE) via density summation.
    //The density of free-surface particles has already been set via density integration and remains untouched!
    for (int lidNodeCol_i = 0; lidNodeCol_i != discret_->NodeColMap()->NumMyElements(); ++lidNodeCol_i)
    {

      ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(PARTICLE_REINITSHIFT==1)
      {
        if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_NONE)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
      }
      else if(PARTICLE_REINITSHIFT==2)
      {
        if(particle_i.freesurfacedata_.color_field_>=1.0)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
      }
      else if(PARTICLE_REINITSHIFT==3)
      {
        particle_i.density_=particle_i.freesurfacedata_.density_sum_;
      }
      else
        dserror("Only the values 1,2 and 3 are possible for PARTICLE_REINITSHIFT!");

      LINALG::Assemble(*density, particle_i.density_, particle_i.gid_, particle_i.owner_);
    }//loop over i
  }

  if(freeSurfaceType==INPAR::PARTICLE::NormalizedReinitialization)
  {
    //loop over particles and re-initialize interior particles (FS_NONE) via density summation
    //and free-surface particles (FS_DIRECT and FS_INDIRECT) via normalized density summation (Shepard Filter)
    for (int lidNodeCol_i = 0; lidNodeCol_i != discret_->NodeColMap()->NumMyElements(); ++lidNodeCol_i)
    {

      ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(PARTICLE_REINITSHIFT==1)
      {
        if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_NONE)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
        else
          particle_i.density_=particle_i.freesurfacedata_.density_sum_/particle_i.freesurfacedata_.color_field_;
      }
      else if(PARTICLE_REINITSHIFT==2)
      {
        if(particle_i.freesurfacedata_.color_field_>=1.0)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
        else
          particle_i.density_=particle_i.freesurfacedata_.density_sum_/particle_i.freesurfacedata_.color_field_;
      }
      else if(PARTICLE_REINITSHIFT==3)
      {
        particle_i.density_=particle_i.freesurfacedata_.density_sum_/particle_i.freesurfacedata_.color_field_;
      }
      else
        dserror("Only the values 1,2 and 3 are possible for PARTICLE_REINITSHIFT!");

      LINALG::Assemble(*density, particle_i.density_, particle_i.gid_, particle_i.owner_);
    }//loop over i
  }

  if(freeSurfaceType==INPAR::PARTICLE::RandlesReinitialization)
  {
    //loop over particles and re-initialize interior particles (FS_NONE) via density summation
    //and free-surface particles (FS_DIRECT and FS_INDIRECT) similar to the procedure in Randles 1996.
    for (int lidNodeCol_i = 0; lidNodeCol_i != discret_->NodeColMap()->NumMyElements(); ++lidNodeCol_i)
    {

      ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(PARTICLE_REINITSHIFT==1)
      {
        if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_NONE)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
        else
        {
          //TODO: Replace p0 with atmospheric pressue in case of inhomogeneous NBC
          double p0=0.0;
          //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
          double density_0=PARTICLE::Utils::Pressure2Density(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),
                  restDensity_,refdensfac_,p0,particle_algorithm_->ExtParticleMat()->exponent_);

          particle_i.density_=particle_i.freesurfacedata_.density_sum_ + density_0*(1.0-particle_i.freesurfacedata_.color_field_);
        }
      }
      else if(PARTICLE_REINITSHIFT==2)
      {
        if(particle_i.freesurfacedata_.color_field_>=1.0)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
        else
        {
          //TODO: Replace p0 with atmospheric pressue in case of inhomogeneous NBC
          double p0=0.0;
          //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
          double density_0=PARTICLE::Utils::Pressure2Density(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),
                  restDensity_,refdensfac_,p0,particle_algorithm_->ExtParticleMat()->exponent_);

          particle_i.density_=particle_i.freesurfacedata_.density_sum_ + density_0*(1.0-particle_i.freesurfacedata_.color_field_);
        }
      }
      else if(PARTICLE_REINITSHIFT==3)
      {
        //TODO: Replace p0 with atmospheric pressue in case of inhomogeneous NBC
        double p0=0.0;
        //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
        double density_0=PARTICLE::Utils::Pressure2Density(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),
                restDensity_,refdensfac_,p0,particle_algorithm_->ExtParticleMat()->exponent_);

        particle_i.density_=particle_i.freesurfacedata_.density_sum_ + density_0*(1.0-particle_i.freesurfacedata_.color_field_);
      }
      else
        dserror("Only the values 1,2 and 3 are possible for PARTICLE_REINITSHIFT!");

      LINALG::Assemble(*density, particle_i.density_, particle_i.gid_, particle_i.owner_);
    }//loop over i
  }

  return;
}

/*--------------------------------------------------------------------------*
 | Initialize free-surface particle - mesh free style  meier         08/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitFreeSurfaceParticles(
    const Teuchos::RCP<Epetra_Vector> fspType)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleMeshFreeInteractionHandler::InitFreeSurfaceParticles");

  freesurfaceparticles_.clear();

  //In the following, particles with color_field < PARTICLE_COLORLIMIT are classified as free-surface particles of type FS_DIRECT. They are summed up
  //in specific data containers freesurfaceparticles_ (containing the particles) and neighbours_fsp_ (containing their neighbors). All neighbors of
  //these particles are classified as free-surface particles of type FS_INDIRECT. These get no extra entries in freesurfaceparticles_ and neighbours_fsp_ so far!

  const int numcolelements = discret_->NodeColMap()->NumMyElements();
  // loop over the particles (no superpositions)
  for (unsigned int lidNodeCol_i = 0; lidNodeCol_i != (unsigned int)(numcolelements); ++lidNodeCol_i)
  {
    // determine the particle_i
    ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    if(particle_i.freesurfacedata_.color_field_ < PARTICLE_COLORLIMIT and particle_i.boundarydata_.boundaryparticle_==false)
    {
      particle_i.freesurfacedata_.freesurfaceparticletype_=FS_DIRECT;
      freesurfaceparticles_[particle_i.gid_]=&colParticles_[lidNodeCol_i];
    }
  }

  //************loop over free surface particles and search for neighbors (for free-surface particles we need superposition --> column map)********
  // bin checker
  std::set<int> examinedbins;

  // loop over the boundary particles
  for (std::map<int,ParticleMF* >::iterator ii = freesurfaceparticles_.begin(); ii!=freesurfaceparticles_.end(); ++ii)
  {
    int id_i = ii->first;
    // extract the node underlying the particle
    DRT::Node *currparticle = discret_->gNode(id_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins2_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // list of particles in Bin
    std::list<DRT::Node*> neighboursLinf_p;

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_DIRECT)
      {
        AddNewNeighbours_fsp(particle_i, neighboursLinf_p, true);
      }
    }
  }
  //***********************************************************************************************************************************************

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeCol_i = 0; lidNodeCol_i != (unsigned int)(numcolelements); ++lidNodeCol_i)
  {
    // determine the particle_i
    ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_INDIRECT)
    {
      if(particle_i.boundarydata_.boundaryparticle_==true)
        dserror("A boundary particle should not be a free-surface particle!");

      freesurfaceparticles_[particle_i.gid_]=&colParticles_[lidNodeCol_i];
    }
  }
  //***********************************************************************************************************************************************

  //************loop over free surface particles and search for neighbors (for free-surface particles we need superposition --> column map)********
  // clear bin checker
  examinedbins.clear();

  // loop over the boundary particles
  for (std::map<int,ParticleMF* >::iterator ii = freesurfaceparticles_.begin(); ii!=freesurfaceparticles_.end(); ++ii)
  {
    int id_i = ii->first;
    // extract the node underlying the particle
    DRT::Node *currparticle = discret_->gNode(id_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins2_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // list of particles in Bin
    std::list<DRT::Node*> neighboursLinf_p;

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_INDIRECT)
      {
        AddNewNeighbours_fsp(particle_i, neighboursLinf_p);
      }
    }
  }
  //***********************************************************************************************************************************************

  // write free-surface particle type in global vector (required for paraview output)
  if(fspType!=Teuchos::null)
  {
    fspType->PutScalar(0.0);
    for (unsigned int lidNodeCol_i = 0; lidNodeCol_i != (unsigned int)(numcolelements); ++lidNodeCol_i)
    {
      // determine the particle_i
      ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      double fspType_i=particle_i.freesurfacedata_.freesurfaceparticletype_;
      LINALG::Assemble(*fspType, fspType_i, particle_i.gid_, particle_i.owner_);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate surface tension interaction potential          sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SurfTensionInterPot(
    const double& radius,
    const double& rRelNorm2,
    double& lambda,
    double& potential
    )
{
  // Note: Following Tartakovsky et al. 2016 the pairwise interaction force F_int_4 is based on a QuinticBspline weight function
  // and the parameter lambda follows from an integration of F_int_4 (see equations (45) and (46) for lambda in three and two dimensions)
  // here we use the normalized QuinticBspline (in contrast to Tartakovsky et al. 2016 and in accordance with Kordilla et al. 2013), hence
  // lambda contains also normalization terms

  Teuchos::RCP<PARTICLE::WeightFunction_QuinticBspline> weightFunction = Teuchos::rcp(new PARTICLE::WeightFunction_QuinticBspline(WF_DIM_));

  // scale factor of repulsive part
  double A = DRT::Problem::Instance()->ParticleParams().get<double>("SURFTENSION_POT_A");
  // scale factor of attractive part
  double B = DRT::Problem::Instance()->ParticleParams().get<double>("SURFTENSION_POT_A");

  // ratio of repulsive and attractive parts
  double ratio = DRT::Problem::Instance()->ParticleParams().get<double>("SURFTENSION_POT_RATIO");

  // range of repulsive part
  double radius0 = ratio*radius;

  switch (WF_DIM_)
  {
  case INPAR::PARTICLE::WF_3D :
  {
    lambda = (7.0*81.0)/(324.0*359.0) * ( -A*std::pow(radius0, 2) + B*std::pow(radius, 2) );
    break;
  }
  case INPAR::PARTICLE::WF_2D :
  {
    lambda = (2771.0*63.0)/(20412.0*478.0*PI) * ( -A*std::pow(radius0, 2) + B*std::pow(radius, 2) );
    break;
  }
  default :
  {
    dserror("Pairwise interaction force currently just possible for dimension 2 and 3!");
  }
  }

  if(lambda <= 0.0)
    dserror("The parameter lambda in the pairwise interaction force should be positive!");

  potential = -A*weightFunction->W(rRelNorm2, radius0) + B*weightFunction->W(rRelNorm2, radius);

  return;
}

/*------------------------------------------------------------------------*
 | return surface tension time fac                           sfuchs 10/17 |
 *------------------------------------------------------------------------*/
double PARTICLE::ParticleMeshFreeInteractionHandler::SurfTensionTimeFac(const double& time)
{
  double fac=1.0;
  double ramp_time=DRT::Problem::Instance()->ParticleParams().get<double>("SURFTENSION_RAMP_TIME");
  if(ramp_time>0 and time >=0)
  {
    if(time<ramp_time)
       fac=0.5*(1-cos(time*PI/ramp_time));
  }

  return fac;
}

/*------------------------------------------------------------------------------*
 | wrapper                                                         katta 01/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::PrintNeighbours(const PARTICLE::NeighbourType& nt)
{
  switch (nt)
  {
  case PARTICLE::NT_Particle :
  {
    PrintNeighbours(neighbours_p_);
    break;
  }
  case PARTICLE::NT_HeatSource :
  {
    PrintNeighbours(neighbours_hs_);
    break;
  }
  }
}
