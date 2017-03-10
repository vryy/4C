/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.cpp

\brief Particle-MeshFree interaction handling.
papers: - Smoothed dissipative particle dynamics, DOI: 10.1103/PhysRevE.67.026705
        - Numerical simulation of fluid-structure interaction by SPH, DOI: 10.1016/j.compstruc.2007.01.002


\level 3

\maintainer Alessandro Cattabiani
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particleMeshFree_interaction.H"
#include "particleMeshFree_weightFunction.H"
#include "particleMeshFree_surfaceTensionInteractions.H"
#include "particle_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "particle_heatSource.H"
#include "particle_utils.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"


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
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE")),
  WF_DIM_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunctionDim>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION_DIM")),
  initDensity_(initDensity),
  restDensity_(restDensity),
  refdensfac_(refdensfac),
  alphaMin_(DRT::Problem::Instance()->ParticleParams().get<double>("ALPHA_MIN")),
  trg_updatedColorFieldGradient_(false)
{
// checks
if (particle_algorithm_->ExtParticleMat() == NULL)
  dserror("extParticleMat_ is empty");
// extract wall parameters
const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();

switch(wallInteractionType_)
{
case INPAR::PARTICLE::InitParticle :
  {
    wallMeshFreeData_.density_ = extParticleMat->initDensity_;
    wallMeshFreeData_.mass_ = extParticleMat->initDensity_ * PARTICLE::Utils::Radius2Volume(extParticleMat->initRadius_);
    // the pressure is linked to the deltaDensity_ with the initial density. In case of this wall it is always 0
    wallMeshFreeData_.pressure_ = 0;
    /*
    if (extParticleMat->initTemperature_ < extParticleMat->transitionTemperature_)
    {
      wallMeshFreeData_.pressure_ = PARTICLE::Utils::Density2Pressure(extParticleMat->SpeedOfSoundS(),extParticleMat->initDensity_);
    }
    else if (extParticleMat->initTemperature_ > extParticleMat->transitionTemperature_)
    {
      wallMeshFreeData_.pressure_ = PARTICLE::Utils::Density2Pressure(extParticleMat->SpeedOfSoundL(),extParticleMat->initDensity_);
    }
    else
    {
      dserror("Start from the transition state not implemented");
    }
    */
    break;
  }
case INPAR::PARTICLE::Mirror :
  {
    wallMeshFreeData_.density_ = -1;
    wallMeshFreeData_.mass_ = -1;
    wallMeshFreeData_.pressure_ = -1;
    break;
  }
case INPAR::PARTICLE::Custom :
  {
    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    wallMeshFreeData_.density_ = particleparams.get<double>("WALL_FAKE_DENSITY");
    wallMeshFreeData_.mass_ = particleparams.get<double>("WALL_FAKE_MASS");
    wallMeshFreeData_.pressure_ = particleparams.get<double>("WALL_FAKE_PRESSURE");
    break;
  }
case INPAR::PARTICLE::BoundarParticle_NoSlip :
case INPAR::PARTICLE::BoundarParticle_FreeSlip :
case INPAR::PARTICLE::NoWallInteraction :
  {
    //nothing to do
  }
}

// other checks
if (wallInteractionType_ == INPAR::PARTICLE::InitParticle or wallInteractionType_ == INPAR::PARTICLE::Custom)
{
  if (wallMeshFreeData_.density_ < 0)
    dserror("the value of WALL_FAKE_DENSITY is unacceptable");
  if (wallMeshFreeData_.mass_ < 0)
    dserror("the value of WALL_FAKE_MASS is unacceptable");
}

//Attention: The trace of a tensor as required to derive equation (28) from equation (25) in Espanol2003 depends
//on the spatial dimension --> SPH approximations of the Laplace operator typically differ for different dimensions!
switch (WF_DIM_)
{
  case INPAR::PARTICLE::WF_3D :
  {
    diffusionCoeff_ = 5.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
    break;
  }
  case INPAR::PARTICLE::WF_2D :
  {
    diffusionCoeff_ = 8.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
    break;
  }
  case INPAR::PARTICLE::WF_1D :
  {
    diffusionCoeff_ = 11.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
    break;
  }
  default :
  {
    dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }
}

convectionCoeff_ = 5.0 * (extParticleMat->dynamicViscosity_ / 3.0  + extParticleMat->bulkViscosity_);

// checks
if (diffusionCoeff_<0)
{
  dserror("The diffusion coefficient is negative! The following equation should hold: 5*dynamicViscosity >= 3*bulkViscosity");
}

if (convectionCoeff_<0)
{
  dserror("The convection coefficient is negative! Are you sure that the dynamic viscosity and the bulk modulus are positive?");
}

surfaceVoidTension_ = extParticleMat->surfaceVoidTension_;
surfaceWallTension_ = extParticleMat->surfaceWallTension_;

// set the correct WeightFunction

switch (DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION"))
{
  case INPAR::PARTICLE::CubicBspline :
  {
    weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_CubicBspline(WF_DIM_));
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

if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM")==false and extParticleMat->initTemperature_<extParticleMat->transitionTemperature_)
  dserror("Pure mechanical problems, i.e. SOLVE_THERMAL_PROBLEM==No, are currently only considered as pure fluid problems. "
      "To remain consistent, the initial temperature has to be higher than the transition temperature!");

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
    Teuchos::RCP<Epetra_Vector> bpDoFs)
{
  // check
  if (colParticles_.size() != 0)
  {
    dserror("you did not call Clear before Init, colParticles_ is not empty");
  }

  // security block
  trg_updatedColorFieldGradient_ = false;

  // set up the local data storage and fill it with the state vectors
  if (!neighbours_p_.empty() || !neighbours_w_.empty() || !neighbours_hs_.empty())
  {
    std::cout << "The neighbours have memory (They are not empty)!\n";
    std::cout << "However, lid row/col ids were not updated (because not yet updated). It is safe only when not parallel\n";
    std::cin.get();
  }
  InitColParticles(bpDoFs);

  // set up positions and radii to set up the neighbours
  SetStateVector(disn, PARTICLE::Dis);
  SetStateVector(radiusn, PARTICLE::Radius);

  // keep going with the remaining state vectors
  SetStateVector(veln, PARTICLE::Vel);
  SetStateVector(mass, PARTICLE::Mass);
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
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<Epetra_Vector> bpDoFs)
{
  Init(step, disn, veln, radiusn, mass, specEnthalpyn,bpDoFs);

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
    Teuchos::RCP<const Epetra_Vector> pressure,
    Teuchos::RCP<Epetra_Vector> bpDoFs)
{
  Init(step, disn, veln, radiusn, mass, specEnthalpyn, temperature, bpDoFs);

  // set the other state vectors
  SetStateVector(densityn, PARTICLE::Density);
  SetStateVector(pressure, PARTICLE::Pressure);
}


/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitColParticles(Teuchos::RCP<Epetra_Vector> bpDoFs)
{
  // row to col vectors
    // dof-based vectors
  //Teuchos::RCP<Epetra_Vector> disnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
  //Teuchos::RCP<Epetra_Vector> velnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
    // node-based vectors
  //Teuchos::RCP<Epetra_Vector> radiusnCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> massCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> densitynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> specEnthalpynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> pressureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> temperatureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  //Teuchos::RCP<Epetra_Vector> densityapproxCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
    // exports
  //LINALG::Export(*disn ,*disnCol);
  //LINALG::Export(*veln ,*velnCol);
  //LINALG::Export(*radiusn ,*radiusnCol);
  //LINALG::Export(*mass ,*massCol);
  //LINALG::Export(*densityn ,*densitynCol);
  //LINALG::Export(*specEnthalpyn ,*specEnthalpynCol);
  //LINALG::Export(*pressure ,*pressureCol);
  //LINALG::Export(*temperature ,*temperatureCol);
  //LINALG::Export(*densityapproxn ,*densityapproxCol);



  const int numcolelements = discret_->NodeColMap()->NumMyElements();

  colParticles_.resize(numcolelements);
  boundaryparticles_.clear();

  Epetra_Vector bpDoFs_col(*(discret_->DofColMap()),true);
  if(bpDoFs!=Teuchos::null)
  {
    // export bpDoFs (which is in row map format) into overlapping column map format, since the boundary particle information
    // of ghosted particles are also required!
    LINALG::Export(*bpDoFs,bpDoFs_col);
  }

  for (int lidNodeCol=0; lidNodeCol<numcolelements; ++lidNodeCol)
  {
    DRT::Node *particle = discret_->lColNode(lidNodeCol);
    std::vector<int> lm;
    lm.reserve(3);
    discret_->Dof(particle, lm);

    bool boundaryparticle = false;
    if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    {
      if(bpDoFs!=Teuchos::null)
      {
        if((bpDoFs_col)[discret_->DofColMap()->LID(lm[0])]==1)
        {
          boundaryparticle = true;
          if((bpDoFs_col)[discret_->DofColMap()->LID(lm[1])]!=1 or (bpDoFs_col)[discret_->DofColMap()->LID(lm[2])]!=1)
            dserror("For boundary particles all three DoFs have to be prescribed by a Dirichlet Condition!");
        }
      }
    }

    // dof-based vectors
    //LINALG::Matrix<3,1> dis, vel;
    //DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disnCol, dis, lm);
    //DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velnCol, vel, lm);
    // node-based vectors
    //const double radius = (*radiusnCol)[lidNodeCol];
    //const double mass = (*massCol)[lidNodeCol];
    //const double density = (*densitynCol)[lidNodeCol];
    //const double specEnthalpy = (*specEnthalpynCol)[lidNodeCol];
    //const double pressure = (*pressureCol)[lidNodeCol];
    //const double temperature = (*temperatureCol)[lidNodeCol];
    //const double densityapprox = (*densityapproxCol)[lidNodeCol];
    // set up the particles


    colParticles_[lidNodeCol] = ParticleMF(particle->Id(), particle->Owner(), lm, boundaryparticle);
    if(boundaryparticle)
    {
      boundaryparticles_[particle->Id()]=&colParticles_[lidNodeCol];
    }
  }

}

/*-------------------------------------------------------------------------------------------------*
 | accelerations, modified pressures, velocities and energy for boundary particles    meier 02/17  |
 *-------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitBoudaryData(Teuchos::RCP<const Epetra_Vector> accn, double &bpintergy)
{



  //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    dserror("Boundary particles are only considered in the context of pure fluid problems so far!");

  // export acceleration accn (which is in row map format) into overlapping column map format, since accelerations of
  // ghosted particles are also required in order to determine averaged accelerations of boundary particles!
  Epetra_Vector accn_col(*(discret_->DofColMap()),true);
  LINALG::Export(*accn,accn_col);

  //Clear energy in the beginning
  bpintergy=0.0;

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;

    ParticleMF *particle_i = boundaryparticles_[id_i];

    //For boundary particles, density and pressure are set to the initial values
    //Some modified definitions of pressure and density required for the calculation of the pressure gradient are calculated below
    particle_i->density_=initDensity_;
    double deltadensity=initDensity_-refdensfac_*restDensity_;

    //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
    particle_i->pressure_=PARTICLE::Utils::Density2Pressure(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),deltadensity);

    //node ID equals particle ID id_i
    DRT::Node *node = discret_->gNode(id_i);
    std::vector<int> dofnode = discret_->Dof(node);
    for(int i=0;i<3;i++)
    {
      particle_i->boundarydata_.acc_(i)=accn_col[discret_->DofColMap()->LID(dofnode[i])];
    }

    LINALG::Matrix<3,1> sumjVjWij(true);
    double sumjWij(0.0);

    double sumjpjWij(0.0);
    LINALG::Matrix<3,1> sumjrhojrijWij(true);

    //TODO: So far, the gravity is assumed to be zero. In case gravity or any other type of body forces are applied,
    // these terms are required in Ref. Adami2012, Eq.(23) and have to be considered in the vector "gravity"!
    LINALG::Matrix<3,1> gravity(true);
    LINALG::Matrix<3,1> gminusacci(true);
    gminusacci.Update(1.0,gravity,1.0);
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
      //sum up the left sum of numerator of Ref. Adami2012, Eq.(22)
      sumjpjWij+=particle_j.pressure_*interData_ij.w_ij_;
      //sum up the right sum of numerator of Ref. Adami2012, Eq.(22)
      sumjrhojrijWij.Update(interData_ij.rRelNorm2_*particle_j.density_*interData_ij.w_ij_,interData_ij.rRelVersor_ij_,1.0);
    }

    //Only set valuces in case there are interacting (=close enough) fluid particles
    if(sumjWij>0)
    {
      //Update boundary particle velocity required for viscous forces, see Ref. Adami2012, Eq.(23)
      LINALG::Matrix<3,1> tileVi(true);
      tileVi.Update(1.0/sumjWij,sumjVjWij,0.0);
      LINALG::Matrix<3,1> vw(true);
      vw.Update(2.0,particle_i->vel_,1.0);

      vw.Update(-1.0,tileVi,1.0);
      particle_i->boundarydata_.velMod_.Update(1.0,vw,0.0);
      //Update boundary particle pressure required for pressure forces, see Ref. Adami2012, Eq.(23)
      double dotProduct = gminusacci.Dot(sumjrhojrijWij);
      double pi=(sumjpjWij+dotProduct)/sumjWij;
      particle_i->boundarydata_.pressureMod_=pi;

      //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
      double density_i=PARTICLE::Utils::Pressure2Density(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),restDensity_,refdensfac_,pi);
      particle_i->boundarydata_.densityMod_=density_i;

      //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
      bpintergy+=PARTICLE::Utils::Density2Energy(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(), density_i, restDensity_, refdensfac_, particle_i->mass_);
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
  // erase neighbours keep memory
  neighbours_p_.clear();
  neighbours_w_.clear();
  neighbours_hs_.clear();
  neighbours_bp_.clear();
  boundaryparticles_.clear();
}


/*--------------------------------------------------------------------------------------*
 | clear data, keep the memory and keep interactions up to step-memory     katta 12/16  |
 *--------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Clear(const int step, const int memory)
{
  // erase colParticles_. keep the memory
  colParticles_.clear();
  // erase neighbours keep memory
  neighbours_hs_.clear();
  neighbours_bp_.clear();
  boundaryparticles_.clear();

  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      if (jj->second.step_ + memory <= step)
      {
        jj = neighbours_p_[lidNodeRow_i].erase(jj);
      }
    }

    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      if (jj->second.step_ + memory <= step)
      {
        jj = neighbours_w_[lidNodeRow_i].erase(jj);
      }
    }
  }

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

/// miraculous transformation into column vector... ///
Teuchos::RCP<Epetra_Vector> stateVectorCol;
switch (svt)
{
// dof based vectors
case PARTICLE::Dis :
case PARTICLE::Vel :
case PARTICLE::ColorFieldGradient :
case PARTICLE::mGradW :
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
  case PARTICLE::ColorFieldGradient :
  {
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.colorFieldGradient_, data.lm_);
    trg_updatedColorFieldGradient_ = true;
    break;
  }
  case PARTICLE::mGradW :
  {
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.mGradW_, data.lm_);
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
    break;
  }
  case PARTICLE::Temperature :
  {
    data.temperature_ = (*stateVectorCol)[lidNodeCol];
    break;
  }
  case PARTICLE::Alpha :
  {
    data.alpha_ = (*stateVectorCol)[lidNodeCol];
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
 | set colVectors in the local data structs                katta 10/16  |
 *----------------------------------------------------------------------*/
/*
void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVector(Teuchos::RCP<const Epetra_MultiVector> stateVector)
{
// checks
if (stateVector == Teuchos::null)
{
  dserror("the state vector is empty");
}
if (stateVector->NumVectors() != 3)
{
  dserror("the multi vector does not contain 3 vectors");
}


/// miraculous transformation into column vector... ///
Teuchos::RCP<Epetra_MultiVector> stateVectorCol;

stateVectorCol = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofColMap()), 3, true));
LINALG::Export(*stateVector ,*stateVectorCol);

// fill particleData_
for (int lidNodeCol=0; lidNodeCol<discret_->NodeColMap()->NumMyElements(); ++lidNodeCol)
{
  ParticleMF& data = colParticles_[lidNodeCol];
  PARTICLE::Utils::ExtractMyValues(*stateVectorCol, data.mHessW_, data.lm_);
}

}
*/

/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours(const int step)
{
  // resize the vectors
  const int numRowParticles = discret_->NodeRowMap()->NumMyElements();
  neighbours_p_.resize(numRowParticles);
  neighbours_w_.resize(numRowParticles);
  neighbours_hs_.resize(numRowParticles);

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

    // do some checks
    if(neighboursLinf_w.size()!=0)
    {
      if(wallInteractionType_!=INPAR::PARTICLE::Mirror and wallInteractionType_!=INPAR::PARTICLE::Custom and wallInteractionType_!=INPAR::PARTICLE::InitParticle)
        dserror("Wall objects defined but no proper interaction law. Is that your intention?");
    }

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      AddNewNeighbours_p(particle_i, neighboursLinf_p, step);

      if(wallInteractionType_==INPAR::PARTICLE::Mirror or wallInteractionType_==INPAR::PARTICLE::InitParticle or wallInteractionType_==INPAR::PARTICLE::Custom)
        AddNewNeighbours_w(particle_i, neighboursLinf_w, step);

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

}


/*----------------------------------------------------------------------*
 | set the neighbours - particles                          katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_p(
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

    if(particle_j.owner_ != myrank_ || (particle_i.gid_ < particle_j.gid_))
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      rRel_ij.Update(1.0, particle_i.dis_, -1.0, particle_j.dis_); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

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
 //     if (w_ij != 0 || w_ji != 0)
  //    {
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
  //    }
    }
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
      rRel_ij.Update(1.0, particle_i.dis_, -1.0, particle_j.dis_); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

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
 | set the neighbours - walls                              katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_w(
    const ParticleMF& particle_i,
    const boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
    const int step)
{
  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn(Teuchos::null);
  Teuchos::RCP<const Epetra_Vector> wallveln(Teuchos::null);
  if(walldiscret != Teuchos::null)
  {
    walldisn = walldiscret->GetState("walldisnp");
    wallveln = walldiscret->GetState("wallvelnp");
  }

  // evaluate contact with walls first
  std::vector<WallInteractionPoint> surfaces;
  std::vector<WallInteractionPoint> lines;
  std::vector<WallInteractionPoint> nodes;

  std::set<int> unusedIds;

  // check whether there is contact between particle i and neighboring walls
  for(boost::unordered_map<int, DRT::Element*>::const_iterator w=neighboursLinf_w.begin(); w!=neighboursLinf_w.end();  ++w)
  {
    DRT::Element* neighboringwallele = w->second;
    const int numnodes = neighboringwallele->NumNode();
    std::vector<int> lm_wall;
    lm_wall.reserve(numnodes * 3);

    std::vector<int> lmowner;
    std::vector<int> lmstride;

    neighboringwallele->LocationVector(*walldiscret,lm_wall,lmowner,lmstride);

    // nodal displacements
    std::vector<double> nodal_disp(numnodes * 3);
    DRT::UTILS::ExtractMyValues(*walldisn,nodal_disp,lm_wall);

    // get current position of nodes: x = X + u
    std::map<int,LINALG::Matrix<3,1> > nodeCoord;
    DRT::Node** wallnodes = neighboringwallele->Nodes();
    for(int counter=0; counter<numnodes; ++counter)
    {
      static LINALG::Matrix<3,1> currpos;
      const double* X = wallnodes[counter]->X();
      currpos(0) = X[0] + nodal_disp[counter*3+0];
      currpos(1) = X[1] + nodal_disp[counter*3+1];
      currpos(2) = X[2] + nodal_disp[counter*3+2];
      nodeCoord[wallnodes[counter]->Id()] = currpos;
    }

    LINALG::Matrix<3,1> nearestPoint;

    //-------find point on wall element with smallest distance to particle_i-------------------
    GEO::ObjectType objecttype = GEO::nearest3DObjectOnElement(neighboringwallele,nodeCoord,particle_i.dis_,nearestPoint);
    //-----------------------------------------------------------------------------------------

    static LINALG::Matrix<3,1> r_i_wall;
    r_i_wall.Update(1.0, nearestPoint, -1.0, particle_i.dis_);
    const double distance_i_wall = r_i_wall.Norm2();
    const double penetration = distance_i_wall-particle_i.radius_;

    if(penetration <= 0.0)
    {
      // get pointer to the current object type of closest point
      std::vector<WallInteractionPoint> *pointer=0;
      switch(objecttype)
      {
      case GEO::SURFACE_OBJECT:
      {
        pointer = &surfaces;
      }
      break;
      case GEO::LINE_OBJECT:
      {
        pointer = &lines;
      }
      break;
      case GEO::NODE_OBJECT:
      {
        pointer = &nodes;
      }
      break;
      default:
        dserror("unknown object type");
      break;
      }

      // check, whether point has already been detected (e.g. one line element between two surfaces)
      bool insert = true;
      for(size_t i=0; i<(*pointer).size(); ++i)
      {
        static LINALG::Matrix<3,1> distance_vector;
        distance_vector.Update(1.0, nearestPoint, -1.0, (*pointer)[i].point_);
        const double distance = distance_vector.Norm2();
        const double adaptedtol = GEO::TOL7 * particle_i.radius_;

        if (distance < adaptedtol)
        {
          // point has already been detected --> do not insert
          insert = false;
          unusedIds.insert(neighboringwallele->Id());
          break;
        }
      }

      // insert contact point with current surface in corresponding map (surf, line, node)
      if(insert)
      {
        (*pointer).push_back(WallInteractionPoint(neighboringwallele->Id(), lm_wall, lmowner, penetration, nearestPoint, nodeCoord));
      }
    }
    // penetration > 0.0 --> contact impossible
    else
      unusedIds.insert(neighboringwallele->Id());
  }

  // find entries of lines and nodes which are within the penetration volume of the current particle
  // hierarchical: surfaces first
  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // within this radius no other contact point can lie: radius = sqrt(r_i^2 - (r_i-|g|)^2)
    const double rminusg = particle_i.radius_ -std::abs(surfaces[s].penetration_);
    const double radius_surface = sqrt(particle_i.radius_ * particle_i.radius_ - rminusg * rminusg);

    for(size_t l=0; l<lines.size(); ++l)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, surfaces[s].point_, -1.0, lines[l].point_);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_surface)
        unusedIds.insert(lines[l].elemId_);
    }
    for(size_t p=0; p<nodes.size(); ++p)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, surfaces[s].point_, -1.0, nodes[p].point_);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_surface)
        unusedIds.insert(nodes[p].elemId_);
    }
  }
  // find entries of nodes which are within the penetration volume of the current particle
  // hierarchical: lines next
  for(size_t l=0; l<lines.size(); ++l)
  {
    // radius = sqrt(r_i^2 - (r_i-|g|)^2)
    const double rminusg = particle_i.radius_ - std::abs(lines[l].penetration_);
    const double radius_line = sqrt(particle_i.radius_ * particle_i.radius_ - rminusg*rminusg);

    for(size_t p=0; p<nodes.size(); ++p)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, lines[l].point_, -1.0, nodes[p].point_);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_line)
        unusedIds.insert(nodes[p].elemId_);
    }
  }

  // write entries of lines and nodes to surfaces if contact has to be evaluated
  for(size_t l=0; l<lines.size(); ++l)
    if( !unusedIds.count(lines[l].elemId_) )
      surfaces.push_back(lines[l]);
  for(size_t p=0; p<nodes.size(); ++p)
    if( !unusedIds.count(nodes[p].elemId_) )
      surfaces.push_back(nodes[p]);

  // evaluate contact between particle_i and entries of surfaces
  const int lidNodeCol_i = discret_->NodeColMap()->LID(particle_i.gid_);
  std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(discret_->lColNode(lidNodeCol_i))->Get_history_wall();
  if(history_wall.size() > 3)
    dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    WallInteractionPoint wallcontact = surfaces[s];

    // create the data that we have to insert
    LINALG::Matrix<3,1> rRel;
    rRel.Update(1.0, particle_i.dis_, -1.0, wallcontact.point_); // inward vector
    const double rRelNorm2 = rRel.Norm2();
    const double w = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
    const double dw = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
    const double ddw = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

    // insert in the set if appropriate
//    if (dw != 0)
//    {
      LINALG::Matrix<3,1> rRelVersor(rRel);
      rRelVersor.Scale(1/rRelNorm2);

      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
      InterDataPvW& interData_ij = (neighbours_w_[lidNodeRow_i])[wallcontact.elemId_];
      interData_ij = InterDataPvW(
        rRelVersor,
        rRelNorm2,
        step,
        w,
        dw,
        ddw);
//    }
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
 | update weights and distances in all the neighbours                 katta 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights(
    const int step)
{
  UpdateWeights_p(step);
  UpdateWeights_w(step);

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    UpdateWeights_bp(step);
}


/*---------------------------------------------------------------------------------*
 | update weights and distances in all the neighbours - particles     katta 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights_p(const int step)
{
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // shall I skip?
      const int& lidNodeCol_j = jj->first;
      InterDataPvP& interData_ij =(neighbours_p_[lidNodeRow_i])[lidNodeCol_j];
      if (step != interData_ij.step_)
      {
        // extract data for faster access
        const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

        // --- extract and compute general data --- //

        // create the data that we have to insert
        LINALG::Matrix<3,1> rRel_ij;
        rRel_ij.Update(1.0, particle_i.dis_, -1.0, particle_j.dis_); // inward vector
        const double rRelNorm2 = rRel_ij.Norm2();
        LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
        rRelVersor_ij.Scale(1/rRelNorm2);

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


        // insert the new weights, do not touch the step
        interData_ij = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          interData_ij.step_,
          w_ij,
          w_ji,
          dw_ij,
          dw_ji,
          ddw_ij,
          ddw_ji);
      }
    }
  }
}

/*------------------------------------------------------------------------------------------*
 | update weights and distances in all the neighbours - boundary particles     meier 02/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights_bp(const int step)
{

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;
    const ParticleMF *particle_i = boundaryparticles_[id_i];

    // loop over the boundary particles
    for (std::map<int, InterDataPvP>::iterator jj = neighbours_bp_[id_i].begin(); jj!=neighbours_bp_[id_i].end(); ++jj)
    {
      int id_j = jj->first;
      InterDataPvP& interData_ij =(neighbours_bp_[id_i])[id_j];

      if (step != interData_ij.step_)
      {
        // extract data for faster access
        int lidNodeCol_j = discret_->NodeColMap()->LID(id_j);;
        const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

        // --- extract and compute general data --- //

        // create the data that we have to insert
        LINALG::Matrix<3,1> rRel_ij;
        rRel_ij.Update(1.0, particle_i->dis_, -1.0, particle_j.dis_); // inward vector
        const double rRelNorm2 = rRel_ij.Norm2();
        LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
        rRelVersor_ij.Scale(1/rRelNorm2);

        #ifdef PARTICLE_TENSILESAFETYFAC
          //Check that the particle distance does not become to small (danger of tensile instabilities).
          //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
          double initAverageDist = PARTICLE::Utils::Volume2EffDist(std::max(particle_i->mass_,particle_j.mass_)/initDensity_,WF_DIM_);
          if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
            dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
        #endif

        const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i->radius_);
        const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i->radius_);
        const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i->radius_);

        double w_ji = 0;
        double dw_ji = 0;
        double ddw_ji = 0;

        if (particle_i->radius_ == particle_j.radius_)
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

        // insert the new weights, do not touch the step
        interData_ij = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          interData_ij.step_,
          w_ij,
          w_ji,
          dw_ij,
          dw_ji,
          ddw_ij,
          ddw_ji);
      }
    }
  }
}

/*---------------------------------------------------------------------------------*
 | update weights and distances in all the neighbours - particles     katta 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights_w(const int step)
{
  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn(Teuchos::null);
  if(walldiscret != Teuchos::null)
  {
    walldisn = walldiscret->GetState("walldisnp");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_w_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      // shall I skip?
      const int& wallId_j = jj->first;
      InterDataPvW& interData_ij =(neighbours_w_[lidNodeRow_i])[wallId_j];

      if (step != interData_ij.step_)
      {
        // --- extract and compute general data --- //

        DRT::Element* elem_j = walldiscret->gElement(wallId_j);
        const int numnodes = elem_j->NumNode();

        // find point on wall element with smallest distance to particle_i
        std::vector<int> lm_wall;
        lm_wall.reserve(numnodes * 3);
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        elem_j->LocationVector(*walldiscret,lm_wall,lmowner,lmstride);
          // nodal displacements
        std::vector<double> nodal_disp(numnodes * 3);
        DRT::UTILS::ExtractMyValues(*walldisn,nodal_disp,lm_wall);
          // get current position of nodes: x = X + u
        std::map<int,LINALG::Matrix<3,1> > nodeCoord;
        DRT::Node** wallnodes = elem_j->Nodes();
        for(int counter=0; counter<numnodes; ++counter)
        {
         static LINALG::Matrix<3,1> currpos;
         const double* X = wallnodes[counter]->X();
         currpos(0) = X[0] + nodal_disp[counter*3+0];
         currpos(1) = X[1] + nodal_disp[counter*3+1];
         currpos(2) = X[2] + nodal_disp[counter*3+2];
         nodeCoord[wallnodes[counter]->Id()] = currpos;
        }

        LINALG::Matrix<3,1> nearestPoint;
        GEO::nearest3DObjectOnElement(elem_j,nodeCoord,particle_i.dis_,nearestPoint);

        LINALG::Matrix<3,1> rRel;
        rRel.Update(1.0, particle_i.dis_, -1.0, nearestPoint); // inward vector
        const double rRelNorm2 = rRel.Norm2();
        LINALG::Matrix<3,1> rRelVersor(rRel);
        rRelVersor.Scale(1/rRelNorm2);
        const double w = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
        const double dw = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
        const double ddw = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

        // insert the new weights, do not touch the step
        interData_ij = InterDataPvW(
            rRelVersor,
            rRelNorm2,
            interData_ij.step_,
            w,
            dw,
            ddw);
      }
    }
  }
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
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // assemble and write
        LINALG::Matrix<3,1> gradW = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ij_);
        double densityDotn_ij = gradW.Dot(vRel_ij) * particle_j.mass_;
        densityDotn_ij *= extMulti;
        LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid_, particle_i.owner_);
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
    }
  }

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
      vRel_ij.Update(1.0, particle_i->vel_, -1.0, particle_j.vel_);

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
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_acc(
    const Teuchos::RCP<Epetra_Vector> accn,
    const double extMulti,
    const bool clcPressure,
    const double time)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
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

    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;

    // determine some useful quantities
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;


    //If required, add artificial viscous damping force in order to determine static equilibrium configurations
    double damping_factor=DRT::Problem::Instance()->ParticleParams().get<double>("VISCOUS_DAMPING");
    if(damping_factor>0)
    {
      LINALG::Matrix<3,1> dampingforce(true);
      dampingforce.Update(-damping_factor/particle_i.mass_, particle_i.vel_);
      LINALG::Assemble(*accn, dampingforce, particle_i.lm_, particle_i.owner_);
    }

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

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);

      const double rho2 = density_i * density_j;
      const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

      LINALG::Matrix<3,1> momentum_ij;
        // pressure
      // compute the pressure term
      //(based on Monaghan2005, Eq(3.18); this expression differs from the variant in Espanol2003, Eq(22) in the general case m_i \neq m_j.
      //However, based on a similar derivation as in Espanol2003, Eq(18-22) it has been verified that also this variant can guarantee
      //for energy conservation of the dissipation-free, time-continous problem, if Monaghan2005, Eq(2.18) is applied for density integration [and not Eq(2.17)!])
      if (clcPressure)
      {
        const double rhoSquare_j = std::pow(density_j,2);
        const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
          momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);
      }
        // diffusion
      // The following two viscous terms are taken from Espanol2003, Eq(30) and re-expressed in terms of mass densities in an equivalent manner.
      // This is necessary since Espanol2003 only specifies the case m_i=m_j=m.
      const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
      momentum_ij.Update(dC_rho2rRelNorm2, vRel_ij, 1.0);
      // convection
      const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
      momentum_ij.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor_ij, 1.0);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        accn_ij.Scale(extMulti);
        LINALG::Assemble(*accn, accn_ij, particle_i.lm_, particle_i.owner_);
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
      }
    }
  }

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    Inter_bpvp_acc(accn,extMulti,clcPressure);

}


/*----------------------------------------------------------------------*
 | evaluate acceleration - bpvp                            meier 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_bpvp_acc(
    const Teuchos::RCP<Epetra_Vector> accn,
    const double extMulti,
    const bool clcPressure)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
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
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      //The viscous terms require modified velocities of the boundary particles according to Ref. Adami2012, Eq. (23)
      vRel_ij.Update(1.0, particle_i->boundarydata_.velMod_, -1.0, particle_j.vel_);

      LINALG::Matrix<3,1> momentum_ij;
        // pressure
      // compute the pressure term
      //(based on Monaghan2005, Eq(3.8); this expression differs from the variant in Espanol2003, Eq(22) in the general case m_i \neq m_j.
      //However, based on a similar derivation as in Espanol2003, Eq(18-22) it has been verified that also this variant can guarantee
      //for energy conservation of the dissipation-free, time-continous problem, if Monaghan2005, Eq(2.18) is applied for density integration [and not Eq(2.17)!])
      if (clcPressure)
      {

        const double rhoSquare_j = std::pow(density_j,2);
        const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
        momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);
      }

      //Viscous interaction with boundary particles is only required in case of a no-slip boundary condition!
      if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip)
      {
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
      }
    }
  }
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


/*----------------------------------------------------------------------*
 | evaluate color field gradient - pvp                     katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_colorFieldGradient(
    const Teuchos::RCP<Epetra_Vector> colorFieldGradientn,
    const double extMulti)
{
  //checks
  if (colorFieldGradientn == Teuchos::null)
  {
    dserror("colorFieldGradientn is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;
      // --- extract general data --- //
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);

      LINALG::Matrix<3,1> colorFieldGradientGeneral_ij;
      colorFieldGradientGeneral_ij.Update(particle_i.radius_ / density_j,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> colorFieldGradientn_ij;
        colorFieldGradientn_ij.Update(generalCoeff_ij, colorFieldGradientGeneral_ij);
        colorFieldGradientn_ij.Scale(extMulti);
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ij, particle_i.lm_, particle_i.owner_);

      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.dw_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> colorFieldGradientn_ji;
        colorFieldGradientn_ji.Update(generalCoeff_ji, colorFieldGradientGeneral_ij);
        colorFieldGradientn_ji.Scale(-extMulti); // actio = - reactio
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_surfaceTensionCFG(
  const Teuchos::RCP<Epetra_Vector> accn,
  const double extMulti)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }
  if (trg_updatedColorFieldGradient_ == false)
  {
    dserror("the color field gradient was not updated");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];
    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;

    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      //TODO: restDensity=
      const double densityCorrectiveTerm = 2 * restDensity_ / (density_i + density_j);

      // rescaling and assembling
      LINALG::Matrix<3,1> acc_ij;
      acc_ij.Update(1.0, particle_i.colorFieldGradient_, -1.0, particle_j.colorFieldGradient_);
      acc_ij.Scale(- extMulti * densityCorrectiveTerm * surfaceVoidTension_);

      LINALG::Assemble(*accn, acc_ij, particle_i.lm_, particle_i.owner_);

      // rescaling and assembpling
      LINALG::Matrix<3,1> acc_ji;
      acc_ji.Update(-1.0, acc_ij); // actio = - reaction
      LINALG::Assemble(*accn, acc_ji, particle_j.lm_, particle_j.owner_);
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate divergence free pressures - pvp                katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_divFreePressureAcc(
    const Teuchos::RCP<Epetra_Vector> divFreePressureAcc,
    const double dt,
    const double extMulti)
{
  //checks
  if (divFreePressureAcc == Teuchos::null)
  {
    dserror("divFreePressureAcc is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // determine some useful quantities
    const double kv_rho_i = particle_i.densityDot_ * particle_i.alpha_ / (particle_i.density_ * dt);

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      const double kv_rho_j = particle_j.densityDot_ * particle_j.alpha_ / (density_j * dt);

      const double kv_rho = - (kv_rho_i + kv_rho_j);
      LINALG::Matrix<3,1> momentum_ij;
      momentum_ij.Update(kv_rho,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        accn_ij.Scale(extMulti);
        LINALG::Assemble(*divFreePressureAcc, accn_ij, particle_i.lm_, particle_i.owner_);
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
        LINALG::Assemble(*divFreePressureAcc, accn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate constant-density pressures - pvp               katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_constDensityPressureAcc(
    const Teuchos::RCP<Epetra_Vector> constDensityPressureAcc,
    const double dt,
    const double extMulti)
{
  //checks
  if (constDensityPressureAcc == Teuchos::null)
  {
    dserror("constDensityPressureAcc is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // determine some useful quantities
    //TODO: restDensity=
    const double k_rho_i = (particle_i.density_ - restDensity_) * particle_i.alpha_ / (particle_i.density_ * dt * dt);

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];
      const double density_j = particle_j.density_;

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      //TODO: restDensity=
      const double k_rho_j = (density_j - restDensity_) * particle_j.alpha_ / (density_j * dt * dt);

      const double kv_rho = - (k_rho_i + k_rho_j);
      LINALG::Matrix<3,1> momentum_ij;
      momentum_ij.Update(kv_rho,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        accn_ij.Scale(extMulti);
        LINALG::Assemble(*constDensityPressureAcc, accn_ij, particle_i.lm_, particle_i.owner_);
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
        LINALG::Assemble(*constDensityPressureAcc, accn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate acceleration gradient - pvp                    katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_gradAccP(
    const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
    const double &restDensity,
    const double extMulti)
{
  //checks
  if (gradAccP == Teuchos::null)
  {
    dserror("gradAccP is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];
    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;
    const double speedOfSound_i = PARTICLE::Utils::SpeedOfSound(particle_i.specEnthalpy_, particle_algorithm_->ExtParticleMat());

    // useful i variables
    //TODO: restDensity=
    const double gradCoeffP_Rho2_i = std::pow(speedOfSound_i,2) * (2 * restDensity - density_i) / std::pow(density_i,3);
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;
      const double speedOfSound_j = PARTICLE::Utils::SpeedOfSound(particle_j.specEnthalpy_, particle_algorithm_->ExtParticleMat());

      // useful j variables
      //TODO: restDensity=
      const double gradCoeffP_Rho2_j = std::pow(speedOfSound_j,2) * (2 * restDensity - density_j) / std::pow(density_j,3);
      const double rhoSquare_j = std::pow(density_j,2);
      const double p_Rho2_j = particle_j.pressure_/rhoSquare_j;

      // --- gradients and hessians --- //

      LINALG::Matrix<3,1> GradW_ij = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ij_);
      LINALG::Matrix<3,1> mGradW_ij(GradW_ij);
      mGradW_ij.Scale(particle_j.mass_);
      LINALG::Matrix<3,1> mGradW_ji = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ji_);
      mGradW_ji.Scale(-particle_i.mass_); // actio = - reactio
      LINALG::Matrix<3,3> HessW_ij = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
          interData_ij.dw_ij_, interData_ij.rRelNorm2_, interData_ij.ddw_ij_);
      LINALG::Matrix<3,3> mHessW_ij(HessW_ij);
      mHessW_ij.Scale(particle_j.mass_);
      LINALG::Matrix<3,3> mHessW_ji = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
          interData_ij.dw_ji_, interData_ij.rRelNorm2_, interData_ij.ddw_ji_);
      mHessW_ji.Scale(particle_i.mass_); // actio = - reactio ... but twice! (no sign change)

      // assemble --- effect of j on i, main-diagonal block
      if (interData_ij.dw_ij_ != 0)
      {
        // --- diadic products --- //

        // diadic product mGradW_i *  gradRho_i^T
        LINALG::Matrix<3,3> diadicP_mGradWi_gradRhoi;
        diadicP_mGradWi_gradRhoi.MultiplyNT(mGradW_ij, particle_i.mGradW_);

        // diadic product mGradW_i * mGradW_j^T
        LINALG::Matrix<3,3> diadicP_mGradWi_mGradWj;
        diadicP_mGradWi_mGradWj.MultiplyNT(mGradW_ij, mGradW_ji);

        // assemble
        LINALG::Matrix<3,3> gradAccP_iji;
          // 1-term
        gradAccP_iji.Update(- gradCoeffP_Rho2_i, diadicP_mGradWi_gradRhoi, 1.0); // the force is negative!
          // 2-term
        gradAccP_iji.Update( gradCoeffP_Rho2_j,diadicP_mGradWi_mGradWj, 1.0); // the force is negative but we derive with respect of the other variables. This sign needs a check
          // 2-term
        gradAccP_iji.Update(- (p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
        // scale
        gradAccP_iji.Scale(extMulti);
        // convert and write
        std::vector<int> lmowner_i(3, particle_i.owner_);
        Epetra_SerialDenseMatrix Aele(3,3);
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            Aele(dimi,dimj) = gradAccP_iji(dimi,dimj);
          }
        }
        gradAccP->Assemble(-1, Aele, particle_i.lm_, lmowner_i, particle_i.lm_);
      }

      // assemble --- effect of j on i, off-diagonal block
      if (interData_ij.dw_ij_ != 0)
      {
        // --- diadic products --- //

        // diadic product mGradW_i *  gradRho_i^T
        LINALG::Matrix<3,3> diadicP_mGradWi_GradRhoj;
        diadicP_mGradWi_GradRhoj.MultiplyNT(mGradW_ij, particle_j.mGradW_);

        // diadic product mGradW_i * mGradW_j^T
        LINALG::Matrix<3,3> diadicP_mGradWi_mGradWi;
        diadicP_mGradWi_mGradWi.MultiplyNT(mGradW_ij, mGradW_ij);

        // assemble
        LINALG::Matrix<3,3> gradAccP_ijj;
          // 1-term
        gradAccP_ijj.Update(- gradCoeffP_Rho2_j, diadicP_mGradWi_GradRhoj, 1.0); // the force is negative!
          // 2-term
        gradAccP_ijj.Update( gradCoeffP_Rho2_j,diadicP_mGradWi_mGradWi, 1.0); // the force is negative but we derive with respect of the other variables. This sign needs a check
          // 2-term
        gradAccP_ijj.Update( (p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
        // scale
        gradAccP_ijj.Scale(extMulti);
        // convert and write
        std::vector<int> lmowner_i(3, particle_i.owner_);
        Epetra_SerialDenseMatrix Aele(3,3);
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            Aele(dimi,dimj) = gradAccP_ijj(dimi,dimj);
          }
        }
        gradAccP->Assemble(-1, Aele, particle_i.lm_, lmowner_i, particle_j.lm_);
      }

      // assemble --- effect of i on j, main-diagonal block
      if (interData_ij.dw_ji_ != 0 && myrank_ == particle_j.owner_)
      {
        // --- diadic products --- //

        // diadic product mGradW_i *  gradRho_i^T
        LINALG::Matrix<3,3> diadicP_mGradWj_gradRhoj;
        diadicP_mGradWj_gradRhoj.MultiplyNT(mGradW_ji, particle_j.mGradW_);

        // diadic product mGradW_i * mGradW_j^T
        LINALG::Matrix<3,3> diadicP_mGradWj_mGradWi;
        diadicP_mGradWj_mGradWi.MultiplyNT(mGradW_ji, mGradW_ij);

        // assemble
        LINALG::Matrix<3,3> gradAccP_jij;
          // 1-term
        gradAccP_jij.Update(- gradCoeffP_Rho2_i, diadicP_mGradWj_gradRhoj, 1.0); // the force is negative!
          // 2-term
        gradAccP_jij.Update( gradCoeffP_Rho2_j, diadicP_mGradWj_mGradWi, 1.0); // the force is negative but we derive with respect of the other variables. This sign needs a check
          // 2-term
        gradAccP_jij.Update(- (p_Rho2_i+p_Rho2_j), mHessW_ji, 1.0); // the force is negative!
        // scale
        gradAccP_jij.Scale(extMulti);
        // convert and write
        std::vector<int> lmowner_j(3, particle_j.owner_);
        Epetra_SerialDenseMatrix Aele(3,3);
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            Aele(dimi,dimj) = gradAccP_jij(dimi,dimj);
          }
        }
        gradAccP->Assemble(-1, Aele, particle_j.lm_, lmowner_j, particle_j.lm_);
      }

      // assemble --- effect of j on i, off-diagonal block
      if (interData_ij.dw_ji_ != 0 && myrank_ == particle_j.owner_)
      {
        // --- diadic products --- //

        // diadic product mGradW_i *  gradRho_i^T
        LINALG::Matrix<3,3> diadicP_mGradWj_GradRhoi;
        diadicP_mGradWj_GradRhoi.MultiplyNT(mGradW_ji, particle_i.mGradW_);

        // diadic product mGradW_i * mGradW_j^T
        LINALG::Matrix<3,3> diadicP_mGradWj_mGradWj;
        diadicP_mGradWj_mGradWj.MultiplyNT(mGradW_ji, mGradW_ji);

        // assemble
        LINALG::Matrix<3,3> gradAccP_jii;
          // 1-term
        gradAccP_jii.Update(- gradCoeffP_Rho2_j, diadicP_mGradWj_GradRhoi, 1.0); // the force is negative!
          // 2-term
        gradAccP_jii.Update( gradCoeffP_Rho2_j,diadicP_mGradWj_mGradWj, 1.0); // the force is negative but we derive with respect of the other variables. This sign needs a check
          // 3-term
        gradAccP_jii.Update( (p_Rho2_i+p_Rho2_j), mHessW_ji, 1.0); // the force is negative!
        // scale
        gradAccP_jii.Scale(extMulti);
        // convert and write
        std::vector<int> lmowner_j(3, particle_j.owner_);
        Epetra_SerialDenseMatrix Aele(3,3);
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            Aele(dimi,dimj) = gradAccP_jii(dimi,dimj);
          }
        }
        gradAccP->Assemble(-1, Aele, particle_j.lm_, lmowner_j, particle_i.lm_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate acceleration gradient - pvp                    katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_gradAccPapproxOnlyHess(
    const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
    const double &restDensity,
    const double extMulti)
{
  //checks
  if (gradAccP == Teuchos::null)
  {
    dserror("gradAccP is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];
    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;

    // useful i variables
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;


    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      //const double density_j = useDensityApprox ? particle_j.densityapprox_ : particle_j.density_;
      const double density_j = particle_j.density_;

      // useful j variables
      const double rhoSquare_j = std::pow(density_j,2);
      const double p_Rho2_j = particle_j.pressure_/rhoSquare_j;

      // --- gradients and hessians --- //
      LINALG::Matrix<3,3> HessW_ij = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
          interData_ij.dw_ij_, interData_ij.rRelNorm2_, interData_ij.ddw_ij_);
      LINALG::Matrix<3,3> mHessW_ij(HessW_ij);
      mHessW_ij.Scale(particle_j.mass_);
      LINALG::Matrix<3,3> mHessW_ji = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
          interData_ij.dw_ji_, interData_ij.rRelNorm2_, interData_ij.ddw_ji_);
      mHessW_ji.Scale(particle_i.mass_); // actio = - reactio ... but twice! (no sign change)

      // assemble --- effect of j on i, main-diagonal block
      if (interData_ij.dw_ij_ != 0)
      {
        // assemble
        LINALG::Matrix<3,3> gradAccP_iji;
          // 3-term
        gradAccP_iji.Update(- (p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
        // scale
        gradAccP_iji.Scale(extMulti);
        //write
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            gradAccP->Assemble(gradAccP_iji(dimi,dimj), particle_i.lm_[dimi], particle_i.lm_[dimj]);
          }
        }
      }

      // assemble --- effect of j on i, off-diagonal block
      if (interData_ij.dw_ij_ != 0)
      {
        // assemble
        LINALG::Matrix<3,3> gradAccP_ijj;
          // 3-term
        gradAccP_ijj.Update( (p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
        // scale
        gradAccP_ijj.Scale(extMulti);
        // convert and write
        std::vector<int> lmowner_i(3, particle_i.owner_);
        //write
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            gradAccP->Assemble(gradAccP_ijj(dimi,dimj), particle_i.lm_[dimi], particle_j.lm_[dimj]);
          }
        }
      }

      // assemble --- effect of i on j, main-diagonal block
      if (interData_ij.dw_ji_ != 0 && myrank_ == particle_j.owner_)
      {
        // assemble
        LINALG::Matrix<3,3> gradAccP_jij;
          // 3-term
        gradAccP_jij.Update(- (p_Rho2_i+p_Rho2_j), mHessW_ji, 1.0); // the force is negative!
        // scale
        gradAccP_jij.Scale(extMulti);
        //write
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            gradAccP->Assemble(gradAccP_jij(dimi,dimj), particle_j.lm_[dimi], particle_j.lm_[dimj]);
          }
        }
      }

      // assemble --- effect of j on i, off-diagonal block
      if (interData_ij.dw_ji_ != 0 && myrank_ == particle_j.owner_)
      {
        // assemble
        LINALG::Matrix<3,3> gradAccP_jii;
          // 3-term
        gradAccP_jii.Update( (p_Rho2_i+p_Rho2_j), mHessW_ji, 1.0); // the force is negative!
        // scale
        gradAccP_jii.Scale(extMulti);
        //write
        for (int dimi=0; dimi<3; ++dimi)
        {
          for (int dimj=0; dimj<3; ++dimj)
          {
            gradAccP->Assemble(gradAccP_jii(dimi,dimj), particle_j.lm_[dimi], particle_i.lm_[dimj]);
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | evaluate density - pvw                                  katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_densityDot(
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

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const InterDataPvW& interData_ij = jj->second;

      // --- extract general data --- //
      LINALG::Matrix<3,1> rRelVersor(interData_ij.rRelVersor_);
      LINALG::Matrix<3,1> vRel(particle_i.vel_);

      if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
      {
        wallMeshFreeData_.density_ = particle_i.density_;
        wallMeshFreeData_.mass_ = particle_i.mass_;
        // mirrored velocity
        vRel.Scale(2);
      }

      const double rRelVersorDotVrel = rRelVersor.Dot(vRel);
      const double density = rRelVersorDotVrel;

      if (interData_ij.dw_ != 0)
      {
        // construct the specific coeff
        const double generalCoeff = interData_ij.dw_ * wallMeshFreeData_.mass_;
        // assemble and write
        double densityDotn_ij = generalCoeff * density;
        densityDotn_ij *= extMulti;
        LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid_, particle_i.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate acceleration - pvw                             katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_acc(
    const Teuchos::RCP<Epetra_Vector> accn,
    const double extMulti)
{
  //checks
  if (accn == Teuchos::null)
  {
    dserror("accn is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];
    //const double density_i = useDensityApprox ? particle_i.densityapprox_ : particle_i.density_;
    const double density_i = particle_i.density_;

    // determine some useful quantities
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const InterDataPvW& interData_ij = jj->second;

       // --- extract general data --- //

       const double rRelNorm2 = interData_ij.rRelNorm2_;
       LINALG::Matrix<3,1> rRelVersor(interData_ij.rRelVersor_);
       LINALG::Matrix<3,1> vRel(particle_i.vel_);

       if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
       {
         wallMeshFreeData_.pressure_ = particle_i.pressure_;
         wallMeshFreeData_.density_ = density_i;
         wallMeshFreeData_.mass_ = particle_i.mass_;
         // mirrored velocity
         vRel.Scale(2);
       }

       const double rho2 = density_i * wallMeshFreeData_.density_;
       const double rRelVersorDotVrel = rRelVersor.Dot(vRel);

       // --- momentum --- //

       LINALG::Matrix<3,1> momentum;
         // pressure
       const double rhoSquare_j = std::pow(wallMeshFreeData_.density_,2);
       const double gradP_Rho2 = (p_Rho2_i + wallMeshFreeData_.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
       momentum.Update(- gradP_Rho2, rRelVersor, 1.0);

         // diffusion
       const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
       momentum.Update(dC_rho2rRelNorm2, vRel, 1.0);
       // convection
       const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
       momentum.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor, 1.0);

       if (interData_ij.dw_ != 0)
       {
         // construct the specific coeff
         const double generalCoeff = interData_ij.dw_ * wallMeshFreeData_.mass_;
         // assemble and write
         LINALG::Matrix<3,1> accn_ij;
         accn_ij.Update(generalCoeff, momentum);
         accn_ij.Scale(extMulti);
         LINALG::Assemble(*accn, accn_ij, particle_i.lm_, particle_i.owner_);
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

/*----------------------------------------------------------------------*
 | evaluate acceleration gradient - pvp                    katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_gradAccP(
    const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
    const double &restDensity,
    const double extMulti)
{

  //checks
  if (gradAccP == Teuchos::null)
  {
    dserror("gradAccP is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    const double density_i = particle_i.density_;
    const double speedOfSound_i = PARTICLE::Utils::SpeedOfSound(particle_i.specEnthalpy_, particle_algorithm_->ExtParticleMat());

    // this changes with the pressure-density relation
    const double dP_drho_i = speedOfSound_i * speedOfSound_i;

    // useful i variables
    const double rhoSquare_i = std::pow(density_i,2);
    const double coeff_dP_Rho2_drho_i = dP_drho_i / rhoSquare_i - 2 * particle_i.pressure_ / std::pow(density_i,3);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const InterDataPvW& interData_ij = jj->second;

      if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
      {
        wallMeshFreeData_.pressure_ = particle_i.pressure_;
        wallMeshFreeData_.density_ = density_i;
        wallMeshFreeData_.mass_ = particle_i.mass_;
      }

      // useful j variables
      const double rhoSquare_j = std::pow(wallMeshFreeData_.density_,2);
      const double p_Rho2_j = wallMeshFreeData_.pressure_/rhoSquare_j;

       // --- gradients and hessians --- //

       LINALG::Matrix<3,1> GradW_ij = weightFunctionHandler_->GradW(interData_ij.rRelVersor_, interData_ij.dw_);
       LINALG::Matrix<3,1> mGradW_ij(GradW_ij);
       mGradW_ij.Scale(wallMeshFreeData_.mass_);
       LINALG::Matrix<3,3> HessW_ij = weightFunctionHandler_->HessW(interData_ij.rRelVersor_,
           interData_ij.dw_, interData_ij.rRelNorm2_, interData_ij.ddw_);
       LINALG::Matrix<3,3> mHessW_ij(HessW_ij);
       mHessW_ij.Scale(wallMeshFreeData_.mass_);

       // assemble --- effect of j on i, main-diagonal block
       if (interData_ij.dw_ != 0 || interData_ij.ddw_ != 0)
       {
         // --- diadic products --- //

         // diadic product mGradW_i *  gradRho_i^T
         LINALG::Matrix<3,3> diadicP_mGradWi_GradRhoi;
         diadicP_mGradWi_GradRhoi.MultiplyNT(mGradW_ij, particle_i.mGradW_);

         // assemble
         LINALG::Matrix<3,3> gradAccP_iji;
           // 1-term
         double mirrorMulti = 1.0;
         if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
         {
           mirrorMulti = 2.0;
         }
         gradAccP_iji.Update(-mirrorMulti * coeff_dP_Rho2_drho_i, diadicP_mGradWi_GradRhoi, 1.0); // the force is negative!
           // 3-term
         gradAccP_iji.Update(-(p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
         // scale
         gradAccP_iji.Scale(extMulti);
         // convert and write
         std::vector<int> lmowner_i(3, particle_i.owner_);
         Epetra_SerialDenseMatrix Aele(3,3);
         for (int dimi=0; dimi<3; ++dimi)
         {
           for (int dimj=0; dimj<3; ++dimj)
           {
             Aele(dimi,dimj) = gradAccP_iji(dimi,dimj);
           }
         }
         gradAccP->Assemble(-1, Aele, particle_i.lm_, lmowner_i, particle_i.lm_);
       }
    }
  }
}

/*----------------------------------------------------------------------*
 | evaluate acceleration gradient - pvp                    katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_gradAccPapproxOnlyHess(
    const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
    const double &restDensity,
    const double extMulti)
{

  //checks
  if (gradAccP == Teuchos::null)
  {
    dserror("gradAccP is empty");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    const double density_i = particle_i.density_;

    // useful i variables
    const double rhoSquare_i = std::pow(density_i,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin(); jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const InterDataPvW& interData_ij = jj->second;

      if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
      {
        wallMeshFreeData_.pressure_ = particle_i.pressure_;
        wallMeshFreeData_.density_ = density_i;
        wallMeshFreeData_.mass_ = particle_i.mass_;
      }

      // useful j variables
      const double rhoSquare_j = std::pow(wallMeshFreeData_.density_,2);
      const double p_Rho2_j = wallMeshFreeData_.pressure_/rhoSquare_j;

       // --- hessian --- //
       LINALG::Matrix<3,3> HessW_ij = weightFunctionHandler_->HessW(interData_ij.rRelVersor_,
           interData_ij.dw_, interData_ij.rRelNorm2_, interData_ij.ddw_);
       LINALG::Matrix<3,3> mHessW_ij(HessW_ij);
       mHessW_ij.Scale(wallMeshFreeData_.mass_);

       // assemble --- effect of j on i, main-diagonal block
       if (interData_ij.dw_ != 0 || interData_ij.ddw_ != 0)
       {
         // assemble
         LINALG::Matrix<3,3> gradAccP_iji;
           // 3-term-approx!
         gradAccP_iji.Update(-(p_Rho2_i+p_Rho2_j), mHessW_ij, 1.0); // the force is negative!
         // scale
         gradAccP_iji.Scale(extMulti);
         // convert and write
         std::vector<int> lmowner_i(3, particle_i.owner_);
         Epetra_SerialDenseMatrix Aele(3,3);
         for (int dimi=0; dimi<3; ++dimi)
         {
           for (int dimj=0; dimj<3; ++dimj)
           {
             Aele(dimi,dimj) = gradAccP_iji(dimi,dimj);
           }
         }
         gradAccP->Assemble(-1, Aele, particle_i.lm_, lmowner_i, particle_i.lm_);
       }
    }
  }
}


/*--------------------------------------------------------------------------*
 | compute density - mesh free style                           katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_alpha()
{

  //TODO: @katta: What is alpha, where is it required, on which paper / formulas is it based???

  // create the buffer vectors in row mode
  // final alpha vector
  Teuchos::RCP<Epetra_Vector> alpha = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
  // buffer vector to compute | \Sum m grad(W) | ^2
  Teuchos::RCP<Epetra_Vector> alphaDof = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.dw_ij_ != 0)
      {
        // construct the specific ij coeff
         const double generalCoeff_ij = interData_ij.dw_ij_ * particle_j.mass_;
         // assemble and write dof
         LINALG::Matrix<3,1> alphaDof_ij;
         alphaDof_ij.Update(generalCoeff_ij, rRelVersor_ij);
         LINALG::Assemble(*alphaDof, alphaDof_ij, particle_i.lm_, particle_i.owner_);
        // assemble and write node
        double alpha_ij = alphaDof_ij.Dot(alphaDof_ij);
        LINALG::Assemble(*alpha, alpha_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_ != 0)
      {
        // construct the specific ij coeff
         const double generalCoeff_ji = interData_ij.dw_ji_ * particle_i.mass_;
         // assemble and write dof
         LINALG::Matrix<3,1> alphaDof_ji;
         alphaDof_ji.Update(-generalCoeff_ji, rRelVersor_ij); // actio = - reactio
         LINALG::Assemble(*alphaDof, alphaDof_ji, particle_i.lm_, particle_i.owner_);
        // assemble and write node
        double alpha_ji = alphaDof_ji.Dot(alphaDof_ji);
        LINALG::Assemble(*alpha, alpha_ji, particle_i.gid_, particle_i.owner_);
      }
    }
  }

  // compute and insert | \Sum m grad(W) | ^2 into alpha
  for (int ii = 0; ii < alpha->MyLength(); ++ii)
  {
    for (int jj = 0; jj < 3 ; ++jj)
    {
      (*alpha)[ii] += (*alphaDof)[3*ii + jj] * (*alphaDof)[3*ii + jj];

      if ((*alpha)[ii] < alphaMin_)
      {
        (*alpha)[ii] = alphaMin_;
      }
    }
  }

  // update the local ColParticles
  SetStateVector(alpha, PARTICLE::Alpha);
}

/*--------------------------------------------------------------------------*
 | compute \sum m * W (usually the density) - mesh free style  katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW(
    const Teuchos::RCP<Epetra_Vector> mW,
    const bool withWalls,
    const double extMulti)
{
  //checks
  if (mW == Teuchos::null)
  {
    dserror("mW is empty");
  }

  // erase the vector
  mW->PutScalar(0.0);

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

        LINALG::Assemble(*mW, mW_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.w_ji_ != 0)
      {
        // assemble and write
        double mW_ji = interData_ij.w_ji_ * particle_i.mass_;
        mW_ji *= extMulti;
        LINALG::Assemble(*mW, mW_ji, particle_j.gid_, particle_j.owner_);
      }
    }
  }

  // include the walls
  if (withWalls)
  {
    // loop over the walls
    for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_w_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      //Boundary particles are not considered here
      if(particle_i.boundarydata_.boundaryparticle_==true)
        dserror("How did you get here? The combination of boundary particles and wall interaction is not possible yet!");

      // loop over the interaction particle list
      boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin();
      for (; jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const InterDataPvW& interData = jj->second;

        if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
        {
          wallMeshFreeData_.mass_ = particle_i.mass_;
        }

        // write on particle i if appropriate specializing the quantities
        if (interData.w_ != 0)
        {
          // assemble and write
          double mW_i = interData.w_ * wallMeshFreeData_.mass_;
          mW_i *= extMulti;
          LINALG::Assemble(*mW, mW_i, particle_i.gid_, particle_i.owner_);
        }
      }
    }
  }

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    MF_mW_Boundary(mW,extMulti);

}

/*------------------------------------------------------------------------------------------*
 | compute \sum m * W (usually the density) - boundary particle contributions  meier 02/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mW_Boundary(
    const Teuchos::RCP<Epetra_Vector> mW,
    const double extMulti)
{
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

    //Since the density of the boundary particles particle_i should remain unchanged, we have to assemble the initial density
    //value again after all DoFs of the vector mW have been cleared above. This procedure is not required for boundary particle density calculation
    //(or position update) via time integration since vanishing entries for densityDot in the boundary particle DoFs automatically yield an unchanged density.
    LINALG::Assemble(*mW, particle_i->density_, particle_i->gid_, particle_i->owner_);

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
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 | compute pressure difference (for debugging) - mesh free style  Meier 02/17  |
 *-----------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_DeltaP(
    const Teuchos::RCP<Epetra_Vector> PressureDiff)
{

  //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    dserror("Boundary particles are only considered in the context of pure fluid problems so far!");

  //checks
  if (PressureDiff == Teuchos::null)
  {
    dserror("SumGradW or SumVGradW is empty");
  }

  // erase the vector
  PressureDiff->PutScalar(0.0);

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    //Caution: This method intends to calculate the pressure difference with respect to the rest pressure=c_0^2restDensity_
    //Consequently, this is the only spot where we use NO factor refdensfac_ in the calculation of deltadensity!!!
    double deltadensity=0.0;
    double pressure_diff=0.0;

    //Boundary particles are considered separately
    if(particle_i.boundarydata_.boundaryparticle_==true)
    {
      //First, mark boundary particles without fluid neighbors
      if(particle_i.boundarydata_.densityMod_==-1000)
        pressure_diff=-1000;
      else
      {
        deltadensity=particle_i.boundarydata_.densityMod_-restDensity_;

        //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
        pressure_diff=PARTICLE::Utils::Density2Pressure(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),deltadensity);
      }
    }
    else
    {
      deltadensity=particle_i.density_-restDensity_;

      //TODO: So far, boundary particles are only considered in the context of fluid problems (use of SpeedOfSoundL())!
      pressure_diff=PARTICLE::Utils::Density2Pressure(particle_algorithm_->ExtParticleMat()->SpeedOfSoundL(),deltadensity);
    }

    LINALG::Assemble(*PressureDiff, pressure_diff, particle_i.gid_, particle_i.owner_);
  }
}


/*------------------------------------------------------------------------------*
 | compute \sum m * gradW - mesh free style                        katta 01/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mGradW(
    const Teuchos::RCP<Epetra_Vector> mGradW,
    const bool withWalls,
    const double extMulti)
{
  //checks
  if (mGradW == Teuchos::null)
  {
    dserror("mGradW is empty");
  }

  // erase the vector
  mGradW->PutScalar(0.0);

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // write on particle i if appropriate specializing the quantities

      if (interData_ij.dw_ij_)
      {
        // assemble and write
        LINALG::Matrix<3,1> mGradW_ij = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ij_);
        mGradW_ij.Scale(particle_j.mass_);
        mGradW_ij.Scale(extMulti);
        LINALG::Assemble(*mGradW, mGradW_ij, particle_i.lm_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_)
      {
        // assemble and write
        LINALG::Matrix<3,1> mGradW_ji = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ji_);
        mGradW_ji.Scale(- particle_i.mass_); /// actio = - reaction... we are using the rRelVersor_ij instead of rRelVersor_ji
        mGradW_ji.Scale(extMulti);
        LINALG::Assemble(*mGradW, mGradW_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }

  // include the walls
  if (withWalls)
  {
    // loop over the walls
    for (size_t lidNodeRow_i = 0; lidNodeRow_i != neighbours_w_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      // loop over the interaction particle list
      boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin();
      for (; jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const InterDataPvW& interData = jj->second;

        if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
        {
          wallMeshFreeData_.mass_ = particle_i.mass_;
        }

        // write on particle i if appropriate specializing the quantities
        if (interData.dw_)
        {
          // assemble and write
          LINALG::Matrix<3,1> mGradW_i = weightFunctionHandler_->GradW(interData.rRelVersor_, interData.dw_);
          mGradW_i.Scale(wallMeshFreeData_.mass_);
          mGradW_i.Scale(extMulti);
          LINALG::Assemble(*mGradW, mGradW_i, particle_i.lm_, particle_i.owner_);
        }
      }
    }
  }
}


/*------------------------------------------------------------------------------*
 | compute \sum m * hessW - mesh free style                        katta 01/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_mHessW(
    const Teuchos::RCP<Epetra_MultiVector> mHessW,
    const bool withWalls,
    const double extMulti)
{
  //checks
  if (mHessW == Teuchos::null)
  {
    dserror("mHessW is empty");
  }

  // erase the vector
  mHessW->PutScalar(0.0);

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // write on particle i if appropriate specializing the quantities

      if (interData_ij.dw_ij_ || interData_ij.ddw_ij_)
      {
        // assemble and write
        LINALG::Matrix<3,3> mHessW_ij = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
            interData_ij.dw_ij_, interData_ij.rRelNorm2_, interData_ij.ddw_ij_);
        mHessW_ij.Scale(particle_j.mass_);
        mHessW_ij.Scale(extMulti);
        PARTICLE::Utils::Assemble(*mHessW, mHessW_ij, particle_i.lm_, particle_i.owner_, myrank_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.dw_ji_!=0 || interData_ij.ddw_ji_!=0 )
      {
        // assemble and write
        LINALG::Matrix<3,3> mHessW_ji = weightFunctionHandler_->HessW(interData_ij.rRelVersor_ij_,
            interData_ij.dw_ji_, interData_ij.rRelNorm2_, interData_ij.ddw_ji_);
        mHessW_ji.Scale(particle_i.mass_); /// actio = reaction, second derivative
        mHessW_ji.Scale(extMulti);
        PARTICLE::Utils::Assemble(*mHessW, mHessW_ji, particle_j.lm_, particle_j.owner_, myrank_);
      }
    }
  }

  // include the walls
  if (withWalls)
  {
  // loop over the particles (no superpositions)
    for (size_t lidNodeRow_i = 0; lidNodeRow_i != neighbours_w_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      // loop over the interaction particle list
      boost::unordered_map<int, InterDataPvW>::const_iterator jj = neighbours_w_[lidNodeRow_i].begin();
      for (; jj != neighbours_w_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const InterDataPvW& interData = jj->second;

        if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
        {
          wallMeshFreeData_.mass_ = particle_i.mass_;
        }

        // write on particle i if appropriate specializing the quantities
        if (interData.dw_ || interData.ddw_)
        {
          // assemble and write
          LINALG::Matrix<3,3> mHessW_i = weightFunctionHandler_->HessW(interData.rRelVersor_,
              interData.dw_, interData.rRelNorm2_, interData.ddw_);
          mHessW_i.Scale(wallMeshFreeData_.mass_);
          mHessW_i.Scale(extMulti);
          PARTICLE::Utils::Assemble(*mHessW, mHessW_i, particle_i.lm_, particle_i.owner_, myrank_);
        }
      }
    }
  }
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
  case PARTICLE::NT_Wall :
  {
    PrintNeighbours(neighbours_w_);
    break;
  }
  case PARTICLE::NT_HeatSource :
  {
    PrintNeighbours(neighbours_hs_);
    break;
  }
  }
}
