/*----------------------------------------------------------------------*/
/*!
\file particle_sph_interaction.cpp

\brief Particle SPH interaction handling

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_sph_interaction.H"

#include "particle_algorithm.H"
#include "particle_utils.H"
#include "particle_node.H"
#include "particle_timint_kickdrift.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "particle_sph_rendering.H"
#include "particle_sph_weightFunction.H"
#include "particle_sph_eos.H"

/*----------------------------------------------------------------------*
 | constructor for particle SPH interaction                meier 10/16  |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleSPHInteractionHandler::ParticleSPHInteractionHandler(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams,
  bool norender) :
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  myrank_(discret->Comm().MyPID()),
  weightFunctionHandler_(Teuchos::null),
  equationOfStateHandler_(0),
  interactionVariant2_(DRT::INPUT::IntegralValue<int>(particledynparams,"CALC_ACC_VAR2")),
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(particledynparams,"WALL_INTERACTION_TYPE")),
  freeSurfaceType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(particledynparams,"FREE_SURFACE_TYPE")),
  WF_DIM_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunctionDim>(particledynparams,"WEIGHT_FUNCTION_DIM")),
  background_pressure_(particledynparams.get<double>("BACKGROUND_PRESSURE")),
  transport_velocity_(DRT::INPUT::IntegralValue<int>(particledynparams,"TRANSPORT_VELOCITY")),
  no_veldiff_term_(DRT::INPUT::IntegralValue<int>(particledynparams,"NO_VELDIFF_TERM")),
  damping_factor_(particledynparams.get<double>("VISCOUS_DAMPING")),
  xsph_dampfac_(particledynparams.get<double>("XSPH_DAMPFAC")),
  xsph_stiffac_(particledynparams.get<double>("XSPH_STIFFAC")),
  surfaceTensionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionType>(particledynparams,"SURFACE_TENSION_TYPE")),
  surfTens_ff_(particledynparams.get<double>("SURFACE_TENSION_FF")),
  surfTens_fs_(particledynparams.get<double>("SURFACE_TENSION_FS")),
  surfTensPotA_(particledynparams.get<double>("SURFTENSION_POT_A")),
  surfTensPotB_(particledynparams.get<double>("SURFTENSION_POT_B")),
  surfTensPotRatio_(particledynparams.get<double>("SURFTENSION_POT_RATIO")),
  min_pvp_dist_(1.0e10),
  min_pvw_dist_(1.0e10),
  min_pressure_(1.0e10),
  max_pressure_(-1.0e10)
{
  // checks

  if(wallInteractionType_==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType_==INPAR::PARTICLE::BoundarParticle_FreeSlip)
  {
    INPAR::PARTICLE::ExtendedGhosting extendedGhosting = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ExtendedGhosting>(particledynparams,"EXTENDED_GHOSTING");
    if(extendedGhosting != INPAR::PARTICLE::BdryParticleGhosting and extendedGhosting != INPAR::PARTICLE::AddLayerGhosting)
      dserror("Extended ghosting is required if boundary particles are applied!");
  }

  if(freeSurfaceType_!=INPAR::PARTICLE::FreeSurface_None and freeSurfaceType_!=INPAR::PARTICLE::TwoPhase and particle_algorithm_->BinStrategy()->HavePBC())
    dserror("Periodic free surface flows not possible so far!");

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

  // extract vector of pointers to material parameters
  std::vector<const MAT::PAR::ParticleMat*> particlemat = particle_algorithm_->ParticleMat();
  if (particlemat.size() == 0)
    dserror("No particle material defined!");

  // loop over all particle material parameters
  for(unsigned int i=0; i<particlemat.size(); ++i)
  {
    // cast to extended particle material parameters
    const MAT::PAR::ExtParticleMat* extparticlemat = dynamic_cast<const MAT::PAR::ExtParticleMat*>(particlemat[i]);
    if(extparticlemat == NULL)
      dserror("Could not extract material parameters for extended particle material!");
    // append extended particle material parameter to vector
    extParticleMat_.push_back(extparticlemat);
  }

  // set the equation of state for each particle material
  for(unsigned int i=0; i<extParticleMat_.size(); ++i)
  {
    switch (DRT::INPUT::IntegralValue<INPAR::PARTICLE::EquationOfState>(particledynparams,"EQUATION_OF_STATE"))
    {
      case INPAR::PARTICLE::GenTait :
      {
        equationOfStateHandler_.push_back(Teuchos::rcp(new PARTICLE::EquationOfState_GenTait(extParticleMat_[i]->SpeedOfSoundL(), extParticleMat_[i]->refdensfac_, extParticleMat_[i]->exponent_)));
        break;
      }
      case INPAR::PARTICLE::IdealGas :
      {
        equationOfStateHandler_.push_back(Teuchos::rcp(new PARTICLE::EquationOfState_IdealGas(extParticleMat_[i]->SpeedOfSoundL())));
        break;
      }
    }
  }

  // initialize rendering handler
  const INPAR::PARTICLE::RenderingType renderingType = DRT::INPUT::IntegralValue<INPAR::PARTICLE::RenderingType>(particledynparams,"RENDERING");
  if ((renderingType == INPAR::PARTICLE::StandardRendering or renderingType == INPAR::PARTICLE::NormalizedRendering) and norender==false)
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

  if(freeSurfaceType_==INPAR::PARTICLE::TwoPhase and interactionVariant2_==false)
  {
    if(extParticleMat_.size()!=2)
      dserror("Two-phase flow needs two particle material definitions!");

    if(interactionVariant2_==false)
      dserror("Two-phase flow only possible with acceleration terms according to variant 2!");

    if(fabs(extParticleMat_[0]->artificialViscosity_-extParticleMat_[1]->artificialViscosity_)>1.0e-12)
      dserror("Artificial viscosities of both phases have to be identical!");
  }
}

/*----------------------------------------------------------------------*
 | set up internal variables for future computations       meier 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Init(
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> density0)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::Init");

  // check
  if (colParticles_.size() != 0)
  {
    dserror("you did not call Clear before Init, colParticles_ is not empty");
  }

  // set up the local data storage and fill it with the state vectors
  if (!neighbours_p_.empty())
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
  if(density0!=Teuchos::null)
    SetStateVector(density0, PARTICLE::Density0);

  // set up the neighbours
  AddNewNeighbours();

  return;
}

/*----------------------------------------------------------------------*
 | set all the neighbours                                  meier 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::InitColParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::InitColParticles");

  const int numcolelements = discret_->NodeColMap()->NumMyElements();
  int min_voidparticle_id = DRT::Problem::Instance()->ParticleParams().get<int>("MIN_VOIDPARTICLE_ID");

  colParticles_.resize(numcolelements);
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

    int phase_color = 0;
    int material_id = 0;
    // The -1 is needed since particle IDs start with 1 in the input file and with zero in the code
    if(min_voidparticle_id>0 and particle->Id()>=(min_voidparticle_id-1))
    {
      phase_color = 1;
      material_id = 1;

      if(boundaryparticle)
        dserror("Boundary particles should have IDs < min_voidparticle_id!");
    }
    else if(boundaryparticle)
      phase_color=-1;

    colParticles_[lidNodeCol] = ParticleSPH(particle->Id(), particle->Owner(), lm, boundaryparticle, phase_color, extParticleMat_[material_id]->initDensity_, extParticleMat_[material_id]);

    if(boundaryparticle)
    {
      boundaryparticles_[particle->Id()]=&colParticles_[lidNodeCol];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set density & mass at the very beginning                meier 10/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::InitDensityAndMass(
    const double particle_volume,
    const Teuchos::RCP<Epetra_Vector> density,
    const Teuchos::RCP<Epetra_Vector> mass,
    const Teuchos::RCP<Epetra_Vector> disn,
    const Teuchos::RCP<Epetra_Vector> radiusn)
{
  if(density==Teuchos::null)
    dserror("Empty vector density!");

  if(mass==Teuchos::null)
    dserror("Empty vector mass!");

  density->PutScalar(0.0);
  mass->PutScalar(0.0);


  //in case of two-phase flow, the initial density is determined via density summation!
  const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE");
  bool initSummation=false;
  if(freeSurfaceType==INPAR::PARTICLE::TwoPhase and DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"DENSITY_SUMMATION")==true)
    initSummation=true;

  const int num_row_nodes = discret_->NodeRowMap()->NumMyElements();
  for (int lidNodeRow_i = 0; lidNodeRow_i != num_row_nodes; ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    double density_i=particle_i.extParticleMat_->initDensity_;
    double mass_i=particle_volume*density_i;

    LINALG::Assemble(*mass, mass_i, particle_i.gid_, particle_i.owner_);

    if(initSummation==false)
      LINALG::Assemble(*density, density_i, particle_i.gid_, particle_i.owner_);
  }

  if(initSummation==true)
  {
    // set up positions and radii to set up the neighbours
    SetStateVector(disn, PARTICLE::Dis);
    SetStateVector(radiusn, PARTICLE::Radius);
    SetStateVector(mass, PARTICLE::Mass);

    // set up the neighbours
    AddNewNeighbours();

    // calculate initial density via density summation
    MF_mW(density);
  }

  return;
}

/*-------------------------------------------------------------------------------------------------*
 | accelerations, modified pressures, velocities and energy for boundary particles    meier 02/17  |
 *-------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::InitBoundaryData(
        Teuchos::RCP<const Epetra_Vector> accn,
        const LINALG::Matrix<3,1>& g)
{
  SetStateVector(accn, PARTICLE::Acc);

  // loop over the boundary particles
  for (std::map<int,std::map<int, InterDataPvP> >::iterator ii = neighbours_bp_.begin(); ii!=neighbours_bp_.end(); ++ii)
  {
    int id_i = ii->first;

    ParticleSPH *particle_i = boundaryparticles_[id_i];

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
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];
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
      double pressure_i=(sumjpjWij+dotProduct)/sumjWij;
      particle_i->boundarydata_.pressureMod_=pressure_i;
    }
    else
    {
      particle_i->boundarydata_.pressureMod_=-1000; //-1000 means there are no relevant fluid particle neighbours for this boundary particle
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | clear data, keep memory                                 meier 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Clear()
{
  // erase colParticles_. keep the memory
  colParticles_.clear();
  // erase neighbours keep memory
  neighbours_p_.clear();
  neighbours_bp_.clear();
  neighbours_fsp_.clear();
  boundaryparticles_.clear();
  freesurfaceparticles_.clear();

  return;
}

/*----------------------------------------------------------------------*
 | set colVectors in the local data structs                meier 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::SetStateVector(Teuchos::RCP<const Epetra_Vector> stateVector, const PARTICLE::StateVectorType svt)
{
  // checks
  if (stateVector == Teuchos::null)
  {
    dserror("the state vector is empty");
  }

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::SetStateVector");

  /// miraculous transformation into column vector... ///
  Teuchos::RCP<Epetra_Vector> stateVectorCol;
  switch (svt)
  {
    // dof based vectors
    case PARTICLE::Dis :
    case PARTICLE::Vel :
    case PARTICLE::VelConv :
    case PARTICLE::Acc :
    case PARTICLE::CFG :
    case PARTICLE::CFG_FF :
    case PARTICLE::CFG_F1S :
    case PARTICLE::CFG_F2S :
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
    ParticleSPH& data = colParticles_[lidNodeCol];

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
      case PARTICLE::CFG :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.freesurfacedata_.CFG_, data.lm_);
        break;
      }
      case PARTICLE::CFG_FF :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.freesurfacedata_.CFG_FF_, data.lm_);
        break;
      }
      case PARTICLE::CFG_F1S :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.freesurfacedata_.CFG_F1S_, data.lm_);
        break;
      }
      case PARTICLE::CFG_F2S :
      {
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*stateVectorCol, data.freesurfacedata_.CFG_F2S_, data.lm_);
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
      case PARTICLE::Density0 :
      {
        data.density0_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::DensityDot :
      {
        data.densityDot_ = (*stateVectorCol)[lidNodeCol];
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
      case PARTICLE::ColorField :
      {
        data.freesurfacedata_.color_field_ = (*stateVectorCol)[lidNodeCol];
        break;
      }
      case PARTICLE::DistWall :
      {
        data.boundarydata_.distWall_ = (*stateVectorCol)[lidNodeCol];
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

  return;
}

/*----------------------------------------------------------------------*
 | set all the neighbours                                  meier 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours");

  // resize the vectors
  const int numRowParticles = discret_->NodeRowMap()->NumMyElements();
  neighbours_p_.resize(numRowParticles);

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

    // first neighbours_ round
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p, &neighboursLinf_w);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      AddNewNeighbours_p(particle_i, neighboursLinf_p);
    }
  }

  //*************loop over the boundary particles (for boundary particles we need superposition --> column map)**************************************************
  // clear bin checker
  examinedbins.clear();

  // loop over the boundary particles
  for (std::map<int,ParticleSPH* >::iterator ii = boundaryparticles_.begin(); ii!=boundaryparticles_.end(); ++ii)
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
      const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      if(particle_i.boundarydata_.boundaryparticle_==true)
      {
        AddNewNeighbours_bp(particle_i, neighboursLinf_p);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set the neighbours - particles                          meier 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_p(
    ParticleSPH& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_p");

  // self-neighbours not allowed
  // insert the interaction only if meaningful

  // loop over the neighbours particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

    if(particle_j.gid_!=particle_i.gid_)//standard case: interaction ij with i != j
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      LINALG::Matrix<3,1> dis_i=particle_i.dis_;
      LINALG::Matrix<3,1> dis_j=particle_j.dis_;

      // have periodic boundary conditions
      if (particle_algorithm_->BinStrategy()->HavePBC())
        ShiftPeriodicBoundaryPair(dis_i,dis_j,particle_i.radius_,particle_j.radius_);

      rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

      //Only evaluate pairs that are close enough.
      //TODO: In future applications (e.g. implicit time integration), where neighbours are not searched for every new configuration
      //this criterion might be to strict!!!
      if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
        continue;

      #ifdef PARTICLE_TENSILESAFETYFAC
        //Check that the particle distance does not become to small (danger of tensile instabilities).
        //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
        double initAverageDist = PARTICLE::Utils::Volume2EffDist(particle_i.mass_/particle_i.density0_,WF_DIM_);
        if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
        {
          dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
        }

      #endif

      //check for minimal distance
      if(rRelNorm2<min_pvp_dist_)
        min_pvp_dist_=rRelNorm2;

      const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
      const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
      const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

      // push_back
      LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
      rRelVersor_ij.Scale(1/rRelNorm2);
      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
      (neighbours_p_[lidNodeRow_i])[lidNodeCol_j] = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          w_ij,
          dw_ij,
          ddw_ij);

      if(particle_i.boundarydata_.boundaryparticle_==false and particle_j.boundarydata_.boundaryparticle_==false and particle_i.freesurfacedata_.phase_color_!=particle_j.freesurfacedata_.phase_color_)
        particle_i.freesurfacedata_.havephasechange_=true;

      if(particle_j.boundarydata_.boundaryparticle_==true)
      {
        particle_i.boundarydata_.haveboundaryneighbor_=true;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set the neighbours - boundary particles                 meier 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_bp(
    const ParticleSPH& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_bp");

  // self-neighbours not allowed, boundary particles not allowed as neighbours

  // loop over the neighbours_ particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

    //boundary particles not allowed as neighbours
    if(particle_j.boundarydata_.boundaryparticle_==false)
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      LINALG::Matrix<3,1> dis_i=particle_i.dis_;
      LINALG::Matrix<3,1> dis_j=particle_j.dis_;

      // have periodic boundary conditions
      if (particle_algorithm_->BinStrategy()->HavePBC())
        ShiftPeriodicBoundaryPair(dis_i,dis_j,particle_i.radius_,particle_j.radius_);

      rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

      //Only evaluate pairs that are close enough.
      //TODO: In future applications (e.g. implicit time integration), where neighbours are not searched for every new configuration
      //this criterion might be to strict!!!
      if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
        continue;

#ifdef PARTICLE_TENSILESAFETYFAC
      //Check that the particle distance does not become to small (danger of tensile instabilities).
      //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
      double initAverageDist = PARTICLE::Utils::Volume2EffDist(particle_i.mass_/particle_i.density0_,WF_DIM_);
      if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
        dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
#endif

      //check for minimal distance
      if(rRelNorm2<min_pvw_dist_)
        min_pvw_dist_=rRelNorm2;

      const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
      const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
      const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

      // push_back
      LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
      rRelVersor_ij.Scale(1/rRelNorm2);
      (neighbours_bp_[particle_i.gid_])[particle_j.gid_] = InterDataPvP(
          rRelVersor_ij,
          rRelNorm2,
          w_ij,
          dw_ij,
          ddw_ij);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set the neighbours - free-surface particles             meier 08/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_fsp(
    const ParticleSPH& particle_i,
    const std::list<DRT::Node*>& neighboursLinf_p,
    const bool mark_fs_indirect)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::AddNewNeighbours_fsp");

  // self-neighbours not allowed, boundary particles not allowed as neighbours

  // loop over the neighbours_ particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

    if(particle_i.gid_==particle_j.gid_)
      continue;

    // create the data that we have to insert
    LINALG::Matrix<3,1> rRel_ij;
    LINALG::Matrix<3,1> dis_i=particle_i.dis_;
    LINALG::Matrix<3,1> dis_j=particle_j.dis_;

    rRel_ij.Update(1.0, dis_i, -1.0, dis_j); // inward vector
    const double rRelNorm2 = rRel_ij.Norm2();

    //Only evaluate pairs that are close enough.
    //TODO: In future applications (e.g. implicit time integration), where neighbours are not searched for every new configuration
    //this criterion might be to strict!!!
    if(rRelNorm2-std::max(particle_i.radius_,particle_j.radius_)>0)
      continue;

#ifdef PARTICLE_TENSILESAFETYFAC
    //Check that the particle distance does not become to small (danger of tensile instabilities).
    //The quantity initAverageDist represents the average initial distance between  particles considering there mass and initial density.
    double initAverageDist = PARTICLE::Utils::Volume2EffDist(particle_i.mass_/particle_i.density0_,WF_DIM_);
    if(rRelNorm2<PARTICLE_TENSILESAFETYFAC*initAverageDist)
      dserror("Particle distance to small: rRelNorm2=%f, initAverageDist=%f. Danger of tensile instability!",rRelNorm2,initAverageDist);
#endif

    const double w_ij = weightFunctionHandler_->W(rRelNorm2, particle_i.radius_);
    const double dw_ij = weightFunctionHandler_->DW(rRelNorm2, particle_i.radius_);
    const double ddw_ij = weightFunctionHandler_->DDW(rRelNorm2, particle_i.radius_);

    // push_back
    LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
    rRelVersor_ij.Scale(1/rRelNorm2);
    (neighbours_fsp_[particle_i.gid_])[particle_j.gid_] = InterDataPvP(
        rRelVersor_ij,
        rRelNorm2,
        w_ij,
        dw_ij,
        ddw_ij);

    //mark particle_j as indirect free surface particle in case in has not yet been marked as direct free-surface particle
    if(mark_fs_indirect and particle_j.freesurfacedata_.freesurfaceparticletype_!=FS_DIRECT and particle_j.boundarydata_.boundaryparticle_==false)
      particle_j.freesurfacedata_.freesurfaceparticletype_=FS_INDIRECT;
  }

  return;
}

/*---------------------------------------------------------------------------------*
 | apply coordinate shift in case of periodic boundary conditions    sfuchs 10/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::ShiftPeriodicBoundaryPair(
    const LINALG::Matrix<3,1>& dis_1,
    LINALG::Matrix<3,1>& dis_2,
    const double& radius_1,
    const double& radius_2)
{
  // in case both positions 1 and 2 are close to two opposite periodic boundaries, the position dis_2 is shifted correspondingly
  // it is sufficient to only shift dis_2 if the terms evaluated do only depend on the relative distance between two positions

  // maximum particle radius
  double max_radius = std::max(radius_1,radius_2);

  // get bounding box dimensions
  LINALG::Matrix<3,2> xaabb = particle_algorithm_->BinStrategy()->XAABB();

  // loop over all spatial directions
  for(unsigned int dim=0; dim<3; ++dim)
  {
    // apply periodic boundary conditions in current spatial direction
    if( particle_algorithm_->BinStrategy()->HavePBC(dim) )
    {
      // periodic length in current spatial direction
      double pbc_length = particle_algorithm_->BinStrategy()->PBCDelta(dim);

      // particles are not allowed to interact directly AND across the periodic boundary
      // NOTE: in DEM the periodic length should be larger than four times the particle radius
      if( pbc_length < 2.0*max_radius )
        dserror("Particles are not allowed to interact directly AND across the periodic boundary. The minimal periodic length has to be larger than two time the maximum particle radius!");

      if( dis_1(dim) >= (xaabb(dim,1)-max_radius) and dis_2(dim) <= (xaabb(dim,0)+max_radius) )
        dis_2(dim) += pbc_length;
      else if( dis_2(dim) >= (xaabb(dim,1)-max_radius) and dis_1(dim) <= (xaabb(dim,0)+max_radius) )
        dis_2(dim) -= pbc_length;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate density - pvp                                  meier 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_densityDot(
    const Teuchos::RCP<Epetra_Vector> densityDotn)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_densityDot");

  //checks
  if (densityDotn == Teuchos::null)
    dserror("densityDotn is empty");

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    // contributions to boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    // init the sum of density variation of particle i (due to particles j)
    double sumj_densityDotn_ij=0.0;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // determine particle j
      const int& lidNodeCol_j = jj->first;
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

      // extract particle interaction data
      const InterDataPvP& interData_ij = jj->second;

      LINALG::Matrix<3,1> vRel_ij(true);
      if(transport_velocity_==false)
        vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);
      else
        vRel_ij.Update(1.0, particle_i.velConv_, -1.0, particle_j.velConv_);

      LINALG::Matrix<3,1> gradW = weightFunctionHandler_->GradW(interData_ij.rRelVersor_ij_, interData_ij.dw_ij_);
      double densityDotn_ij = 0.0;
      if(freeSurfaceType_!=INPAR::PARTICLE::TwoPhase)
        densityDotn_ij = gradW.Dot(vRel_ij) * particle_j.mass_;
      else
      {
        //For 2-Phase flow we take a slightly modified version according to Adami et al. 2012, Eq (6)
        //As mentioned in this reference, for boundary particles the initial volume is taken.
        double volume_j = 0.0;
        if(particle_j.boundarydata_.boundaryparticle_==true)
          volume_j = particle_j.mass_ / particle_j.density0_;
        else
          volume_j = particle_j.mass_ / particle_j.density_;

        densityDotn_ij = gradW.Dot(vRel_ij) * particle_i.density_ * volume_j;
      }

      sumj_densityDotn_ij += densityDotn_ij;
    } // loop over j

    // assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*densityDotn, sumj_densityDotn_ij, particle_i.gid_, particle_i.owner_);

  } //loop over i

  return;
}

/*----------------------------------------------------------------------*
 | evaluate acceleration - pvp                             meier 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_acc(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> trvl_acc,
    const Teuchos::RCP<Epetra_Vector> acc_A,
    const double time)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_acc");

  //checks
  if (accn == Teuchos::null)
    dserror("accn is empty");

  // Modified particle transport velocities either due to constant background pressure according to Adami2013 or due to XSPH according to Monaghan1989
  // The terms associated with the constant background pressure are similar to the real pressure terms
  // but with the pressures p_i and p_j replaced by the background pressure. All the related variables
  // and terms are marked by a prefix p0, e.g. p0_momentum_ij etc.. All terms related to XSPH are marked by a prefix xsph
  if (trvl_acc != Teuchos::null)
  {
    if(background_pressure_<0.0 and xsph_dampfac_<0.0 and xsph_stiffac_<0.0)
      dserror("Negative values of background pressure, xsph_dampfac and xsph_stifffac! Why is the vector trvl_acc handed in for this case? Set flag TRANSPORT_VELOCITY to no!");

    if(transport_velocity_==false)
      dserror("Vector trvl_acc handed in even though TRANSPORT_VELOCITY is set to no?!?");

    if((xsph_dampfac_<0.0 and xsph_stiffac_>0.0) or (xsph_dampfac_>0.0 and xsph_stiffac_<0.0))
      dserror("Either both or none of the xsph parameters have to be set!");
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    // contributions to boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    // init the sum of accelerations acting on particle i (due to particles j and due to artificial viscous damping)
    LINALG::Matrix<3,1> sumj_accn_ij(true);
    LINALG::Matrix<3,1> sumj_trvl_accn_ij(true);
    LINALG::Matrix<3,1> sumj_accn_ij_A(true);

    // ********** damping (optional) **********
    // If required, add artificial viscous damping force in order to determine static equilibrium configurations
    if(damping_factor_>0)
      sumj_accn_ij.Update(-damping_factor_, particle_i.vel_);

    // loop over the interaction particle list
    boost::unordered_map<int, InterDataPvP>::const_iterator jj;
    for (jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // determine particle j
      const int& lidNodeCol_j = jj->first;
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

      // extract particle interaction data
      const InterDataPvP& interData_ij = jj->second;

      // construct the specific coefficient
      double generalCoeff_ij = Inter_generalCoeff_ij(particle_i, particle_j, interData_ij);

      // ********** pressure **********
      Inter_pressure(&sumj_accn_ij, generalCoeff_ij, particle_i, particle_j, interData_ij);

      // ********** transport velocity (optional) **********
      // background pressure term for modified convection velocities/accelerations
      if(background_pressure_>0)
        Inter_backgroundPressure(&sumj_trvl_accn_ij, generalCoeff_ij, particle_i, particle_j, interData_ij);

      // xsph contributions to transport velocity
      if(xsph_dampfac_>0 or xsph_stiffac_>0)
        Inter_xsph(&sumj_trvl_accn_ij, particle_i, particle_j, interData_ij);

      // additional term div(A) in Navier-Stokes equation due to non-material particle convection
      if(transport_velocity_)
        Inter_transportVelocity_divA(&sumj_accn_ij, &sumj_accn_ij_A, generalCoeff_ij, particle_i, particle_j, interData_ij);

      // ********** viscous forces **********
      Inter_laminarViscosity(&sumj_accn_ij, generalCoeff_ij, particle_i, particle_j, interData_ij);

      // ********** artificial viscous forces **********
      if(particle_i.extParticleMat_->artificialViscosity_>0)
        Inter_artificialViscosity(&sumj_accn_ij, particle_i, particle_j, interData_ij);

      // ********** surface tension **********
      if(surfTens_ff_>0.0 and (surfaceTensionType_==INPAR::PARTICLE::ST_VDW_DIRECT or surfaceTensionType_==INPAR::PARTICLE::ST_VDW_INDIRECT))
        Inter_surfaceTension(&sumj_accn_ij, particle_i, particle_j, interData_ij, time);

    } // loop over j

    // assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*accn, sumj_accn_ij, particle_i.lm_, particle_i.owner_);
    if(background_pressure_>0)
      LINALG::Assemble(*trvl_acc, sumj_trvl_accn_ij, particle_i.lm_, particle_i.owner_);
    if(acc_A!=Teuchos::null)
      LINALG::Assemble(*acc_A, sumj_accn_ij_A, particle_i.lm_, particle_i.owner_);

  } // loop over i

  return;
}

/*----------------------------------------------------------------------*
 | Computes Laplace operator for vector v - pvp            meier 03/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_Laplace_x(
    const ParticleSPH& particle_j,
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

  return;
}

/*-----------------------------------------------------------------------------*
 | Computes Gradient of divergence operator for vector v - pvp    meier 03/17  |
 *-----------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_pvp_GradDiv_x(
    const ParticleSPH& particle_j,
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

  return;
}

/*--------------------------------------------------------------------------*
 | compute \sum m * W (usually the density)                    meier 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::MF_mW(
    const Teuchos::RCP<Epetra_Vector> mW,
    const Teuchos::RCP<Epetra_Vector> unity_vec,
    const Teuchos::RCP<Epetra_Vector> unity_vec_grad,
    const Teuchos::RCP<Epetra_Vector> phaseColor)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::MF_mW");

  //checks
  if (mW == Teuchos::null)
    dserror("mW is empty");

  // erase the vectors
  mW->PutScalar(0.0);
  if(unity_vec!=Teuchos::null)
    unity_vec->PutScalar(0.0);
  if(unity_vec_grad!=Teuchos::null)
    unity_vec_grad->PutScalar(0.0);
  if(phaseColor!=Teuchos::null)
    phaseColor->PutScalar(0.0);

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    if(phaseColor!=Teuchos::null)
    {
      double phase_color = (double)particle_i.freesurfacedata_.phase_color_;
      LINALG::Assemble(*phaseColor, phase_color, particle_i.gid_, particle_i.owner_);
    }

    // contributions to boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    // auto-interaction
    double mW_ii = weightFunctionHandler_->W0(particle_i.radius_) * particle_i.mass_;
    LINALG::Assemble(*mW, mW_ii, particle_i.gid_, particle_i.owner_);

    double unity_vec_aux=0.0;
    if(unity_vec!=Teuchos::null)
    {
      unity_vec_aux=mW_ii/particle_i.density_;
      LINALG::Assemble(*unity_vec, unity_vec_aux, particle_i.gid_, particle_i.owner_);
    }

    // init the sum of weighted density acting on particle i (due to particles j)
    double sumj_mW_ij=0.0;
    double sumj_m_rho_W_ij=0.0;
    LINALG::Matrix<3,1> sumj_m_rho_dW_ij(true);

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // determine particle j
      const int& lidNodeCol_j = jj->first;
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

      // get density of (boundary) particle j
      bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
      double density_j = particle_j.density_;
      if(boundaryParticle_j)
      {
#ifndef PARTICLE_BOUNDARYDENSITY
        density_j = particle_i.density0_;
#else
        // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
        // the equation of state based on material data of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
        double pressure_j = particle_j.boundarydata_.pressureMod_;
        // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
        const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);

        density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
#endif
      }

      // extract particle interaction data
      const InterDataPvP& interData_ij = jj->second;

      //For density summation, we apply the variant miW_ij proposed by Hu et al. 2006 Eq. (15) (with particle_i.mass_ instead of particle_j.mass_)
      double miW_ij = interData_ij.w_ij_ * particle_i.mass_;
      double mW_ij = interData_ij.w_ij_ * particle_j.mass_;
      sumj_mW_ij += miW_ij;
      if(unity_vec!=Teuchos::null)
        sumj_m_rho_W_ij+=mW_ij/density_j;
      if(unity_vec_grad!=Teuchos::null)
        sumj_m_rho_dW_ij.Update(particle_j.mass_/density_j*interData_ij.dw_ij_,interData_ij.rRelVersor_ij_,1.0);
    } // loop over j

    // assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    LINALG::Assemble(*mW, sumj_mW_ij, particle_i.gid_, particle_i.owner_);
    if(unity_vec!=Teuchos::null)
      LINALG::Assemble(*unity_vec, sumj_m_rho_W_ij, particle_i.gid_, particle_i.owner_);
    if(unity_vec_grad!=Teuchos::null)
      LINALG::Assemble(*unity_vec_grad, sumj_m_rho_dW_ij, particle_i.lm_, particle_i.owner_);

  } // loop over i

  // loop over all boundary particles (not only the one with fluid neighbours!) and set density and color field
  for (std::map<int,ParticleSPH* >::iterator ii = boundaryparticles_.begin(); ii!=boundaryparticles_.end(); ++ii)
  {
    ParticleSPH* particle_i = ii->second;

    // Since the density of the boundary particles particle_i should remain unchanged, we have to assemble the initial density
    // value again after all DoFs of the vector mW have been cleared above. This procedure is not required for boundary particle density calculation
    // (or position update) via time integration since vanishing entries for densityDot in the boundary particle DoFs automatically yield an unchanged density.
    double density = particle_i->density0_;
    LINALG::Assemble(*mW, density, particle_i->gid_, particle_i->owner_);
    if(unity_vec!=Teuchos::null)
    {
      double unity = 1.0;
      LINALG::Assemble(*unity_vec, unity, particle_i->gid_, particle_i->owner_);
    }
  }

  return;
}

/*-------------------------------------------------------------------------------*
 | Prescribe normal vector of particles close to the triple point   meier 11/17  |
 *-------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::SetTriplePointNormal(
    const Teuchos::RCP<Epetra_Vector> CFG_F2S)
{
  //Enforce static contact angle by prescribing normal vector of particles close to the triple point as proposed in Breinlinger et al. 2013.

  if(CFG_F2S == Teuchos::null)
    dserror("CFG_F2S empty!");
  else
    CFG_F2S->PutScalar(0.0);

  bool contact_angle_var2=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"CONTACT_ANGLE_VAR2");

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {

    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    //Boundary particles are not considered here
    if(particle_i.boundarydata_.boundaryparticle_==false and particle_i.freesurfacedata_.havephasechange_==true)
    {

      LINALG::Matrix<3,1> normal1(particle_i.freesurfacedata_.CFG_FF_);
      double norm1=normal1.Norm2();

      if(norm1 > 1.0e-10*particle_i.radius_)
        normal1.Scale(1.0/norm1);
      else
        normal1.Scale(0.0);

      // The contact angle force variant 2 does not required to prescribe the normal vector for particles close to the triple point
      if(contact_angle_var2==true)
      {
        LINALG::Assemble(*CFG_F2S, normal1, particle_i.lm_, particle_i.owner_);
      }
      else
      {
        if(particle_i.boundarydata_.haveboundaryneighbor_==false)
          LINALG::Assemble(*CFG_F2S, normal1, particle_i.lm_, particle_i.owner_);
        else
        {

          double initDist=PARTICLE::Utils::Volume2EffDist(particle_i.mass_/particle_i.density0_,WF_DIM_);
          double dwi=particle_i.boundarydata_.distWall_-initDist;

          double fi=1.0;
          double dmax=particle_i.radius_-PARTICLE_ST_TRIPPLEPOINTSMOOTHFAC*initDist;
          if(dwi<0.0)
            fi=0.0;
          else if(dwi<=dmax)
            fi=dwi/dmax;

          if(fi>1.0000000001)
            dserror("fi must be smaller than 1.0!");

          if(fi<-0.0000000001)
            dserror("fi must be larger than 0.0!");

          LINALG::Matrix<3,1> normal_wi(particle_i.freesurfacedata_.CFG_F1S_);
          double normal_wi_norm=normal_wi.Norm2();
          if(normal_wi_norm>1.0e-10*particle_i.radius_)
            normal_wi.Scale(1.0/normal_wi_norm);
          else
            normal_wi.Scale(0.0);

          LINALG::Matrix<3,1> normal_ti(true);
          normal_ti.Update(1.0,normal1,1.0);
          normal_ti.Update(-normal_wi.Dot(normal1),normal_wi,1.0);
          double norm_ti=normal_ti.Norm2();
          if(norm_ti > 1.0e-10*particle_i.radius_)
            normal_ti.Scale(1.0/norm_ti);
          else
            normal_ti.Scale(0.0);

          LINALG::Matrix<3,1> normal_tli(true);
          // convert static contact angle in radians
          double theta_0 = surfTens_fs_*M_PI/180.0;
          normal_tli.Update(sin(theta_0),normal_ti,1.0);
          normal_tli.Update(-cos(theta_0),normal_wi,1.0);

          LINALG::Matrix<3,1> normal2(true);
          normal2.Update(fi,normal1,1.0);
          normal2.Update(1.0-fi,normal_tli,1.0);
          double norm2=normal2.Norm2();

          if(norm2 > 1.0e-10*particle_i.radius_)
            normal2.Scale(1.0/norm2);
          else
            normal2.Scale(0.0);

          LINALG::Assemble(*CFG_F2S, normal2, particle_i.lm_, particle_i.owner_);
        }
      }
    }
  }

}

/*------------------------------------------------------------------------------------------*
 | Calculate a smoothed version of the color field gradient - mesh free style  meier 10/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::MF_SmoothedCFG(
    const Teuchos::RCP<Epetra_Vector> CFG_FF,
    const Teuchos::RCP<Epetra_Vector> CFG_F1S,
    const Teuchos::RCP<Epetra_Vector> CFG_F2S,
    const Teuchos::RCP<Epetra_Vector> distWall)
{
  //checks
  if(CFG_FF == Teuchos::null)
    dserror("CFG_FF empty!");

  if(CFG_F1S == Teuchos::null)
    dserror("CFG_F1S empty!");

  if(CFG_F2S == Teuchos::null)
    dserror("CFG_F2S empty!");

  if(distWall == Teuchos::null)
    dserror("distWall!");

  const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE");
  const INPAR::PARTICLE::SurfaceTensionType surfaceTensionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionType>(DRT::Problem::Instance()->ParticleParams(),"SURFACE_TENSION_TYPE");

  if(freeSurfaceType!=INPAR::PARTICLE::TwoPhase)
  {
    // loop over the particles (no superpositions)
    for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      //Boundary particles are not considered here

      if(particle_i.boundarydata_.boundaryparticle_==true)
        continue;

      LINALG::Matrix<3,1> sumj_m_rho_dW_ij_FF(true);

      // loop over the interaction particle list
      for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const int& lidNodeCol_j = jj->first;
        const InterDataPvP& interData_ij = jj->second;
        const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

        //Boundary particles are not considered here
        if(particle_j.boundarydata_.boundaryparticle_==true)
          continue;

        // see Morris et al. 2000 for the definition of this smoothed color field gradient!
        // write on particle i if appropriate specializing the quantities
        if (interData_ij.w_ij_ != 0)
        {
          if(CFG_FF!=Teuchos::null)
          {
            sumj_m_rho_dW_ij_FF.Update(particle_j.mass_/particle_j.density_*(particle_j.freesurfacedata_.color_field_-particle_i.freesurfacedata_.color_field_)*interData_ij.dw_ij_,interData_ij.rRelVersor_ij_,1.0);
          }
        }
      }//loop over j

      //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
      if(CFG_FF!=Teuchos::null)
        LINALG::Assemble(*CFG_FF, sumj_m_rho_dW_ij_FF, particle_i.lm_, particle_i.owner_);

    }//loop over i
  }
  else if(freeSurfaceType==INPAR::PARTICLE::TwoPhase and surfaceTensionType==INPAR::PARTICLE::ST_CONTI_ADAMI)
  {
    // loop over the particles (no superpositions)
    for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      //Sum up interactions before assembly
      LINALG::Matrix<3,1> sumj_m_rho_dW_ij_FF(true);
      LINALG::Matrix<3,1> sumj_m_rho_dW_ij_F1S(true);
      LINALG::Matrix<3,1> sumj_m_rho_dW_ij_F2S(true);

      // loop over the interaction particle list
      for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const int& lidNodeCol_j = jj->first;
        const InterDataPvP& interData_ij = jj->second;
        const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

        // see Morris et al. 2000 for the definition of this smoothed color field gradient!
        // write on particle i if appropriate specializing the quantities
        if (interData_ij.w_ij_ != 0)
        {
          if(particle_i.boundarydata_.boundaryparticle_==false and particle_j.boundarydata_.boundaryparticle_==false)
          {
            double tilde_c_ji = 0.0;
            if(particle_i.freesurfacedata_.phase_color_!=particle_j.freesurfacedata_.phase_color_ and particle_j.freesurfacedata_.phase_color_==0)
              tilde_c_ji = particle_i.density_/(particle_i.density_+particle_j.density_);

            if(particle_i.freesurfacedata_.phase_color_!=particle_j.freesurfacedata_.phase_color_ and particle_i.freesurfacedata_.phase_color_==0)
              tilde_c_ji = -particle_i.density_/(particle_i.density_+particle_j.density_);

            double V_i=particle_i.mass_/particle_i.density_;
            double V_j=particle_j.mass_/particle_j.density_;

            if(CFG_FF!=Teuchos::null)
            {
              sumj_m_rho_dW_ij_FF.Update((V_i*V_i+V_j*V_j)/V_i*tilde_c_ji*interData_ij.dw_ij_,interData_ij.rRelVersor_ij_,1.0);
            }
          }

          //interaction of (valid) solid particle j with fluid particle i
          if(particle_j.boundarydata_.boundaryparticle_==true and particle_i.boundarydata_.boundaryparticle_==false)
          {
            //solid/boundary particle -> initial volume (identical for all particles)
            double V_j=particle_i.mass_/particle_i.density0_; //V_j=1/sigma_j
            //fluid particle -> current volume
            double V_i=particle_i.mass_/particle_i.density_; //V_i=1/sigma_i

            //For the Adami variant, all boundary particle interaction is assembled into F1S
            sumj_m_rho_dW_ij_F1S.Update(-V_j*V_j/V_i*interData_ij.dw_ij_,interData_ij.rRelVersor_ij_,1.0);
          }
        }
      }//loop over j

    //Assemble all contributions to particle i (we do this at once outside the j-loop since the assemble operation is expensive!)
    if(CFG_FF!=Teuchos::null)
      LINALG::Assemble(*CFG_FF, sumj_m_rho_dW_ij_FF, particle_i.lm_, particle_i.owner_);

    if(particle_i.boundarydata_.haveboundaryneighbor_==true and particle_i.boundarydata_.boundaryparticle_==false)
    {
      LINALG::Matrix<3,1> normalWall=sumj_m_rho_dW_ij_F1S;
      double norm=normalWall.Norm2();
      if(norm>1.0e-10*particle_i.radius_)
        normalWall.Scale(1.0/norm);
      else
        continue;

      if(CFG_F1S!=Teuchos::null)
        LINALG::Assemble(*CFG_F1S, sumj_m_rho_dW_ij_F1S, particle_i.lm_, particle_i.owner_);

      double distWall_aux=1000*particle_i.radius_;

      // loop over the interaction particle list
      for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const int& lidNodeCol_j = jj->first;
        const InterDataPvP& interData_ij = jj->second;
        const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

        if(particle_j.boundarydata_.boundaryparticle_==false)
          continue;

        double dist=interData_ij.rRelNorm2_*normalWall.Dot(interData_ij.rRelVersor_ij_);

        if(dist<distWall_aux)
          distWall_aux=dist;
      }

      if(distWall!=Teuchos::null)
        LINALG::Assemble(*distWall, distWall_aux, particle_i.gid_, particle_i.owner_);
    }

    }//loop over i
  }

  return;
}


/*------------------------------------------------------------------------------------------*
 | Extrapolate CFG to triple point particles - mesh free style                meier 11/17  |
 *------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::ExtrapolateTriplePointCFG(
    const Teuchos::RCP<Epetra_Vector> CFG_FF)
{
  //checks
  if(CFG_FF == Teuchos::null)
    dserror("CFG_FF empty!");
  else
    CFG_FF->PutScalar(0.0);

  //The basic idea of this extrapolation of the CFG to the triple point particles is, to modify the direction of CFG (unit CFG),
  //but to keep the original magnitude.

  double factor_i=2.0; // /in [0 (max smooth), 3 (jump)] 2 is the actual maximum
  double factor_j=2.0; // /in [0 (max smooth), 3 (jump)] 2 is the actual maximum

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    LINALG::Matrix<3,1> CFG_FF_i(particle_i.freesurfacedata_.CFG_FF_);

    //Only particles close to the triple point are modified, the others keep their old CFG_FF value
    if(particle_i.freesurfacedata_.havephasechange_==true and particle_i.boundarydata_.haveboundaryneighbor_==true and particle_i.boundarydata_.boundaryparticle_==false)
    {
      double norm_CFG_FF_i(CFG_FF_i.Norm2());

      // If the norm is very small, the particle is at the fringe of the phase transition region and has negligible influence.
      // Instead of a very small number, nothing is assembled at all.
      if(norm_CFG_FF_i< 1.0e-12/particle_i.radius_)
        continue;

      double initDist=PARTICLE::Utils::Volume2EffDist(particle_i.mass_/particle_i.density0_,WF_DIM_);
      double dwi=particle_i.boundarydata_.distWall_-initDist;
      double fi=1.0;
      double dmaxi=std::max(particle_i.radius_-factor_i*initDist,0.0);
      if(dwi<0.0)
        fi=0.0;
      else if(dwi<=dmaxi)
        fi=dwi/dmaxi;

#ifdef PARTICLE_ST_CFGEXTRAPLATION
      LINALG::Matrix<3,1> sumj_Vj_CFGFFj_W_ij(true);
      double sumj_Vj_W_ij=0.0;

      // loop over the interaction particle list
      for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
      {
        // extract data for faster access
        const int& lidNodeCol_j = jj->first;
        const InterDataPvP& interData_ij = jj->second;
        const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];

        // see Morris et al. 2000 for the definition of this smoothed color field gradient!
        // write on particle i if appropriate specializing the quantities
        if (interData_ij.w_ij_ != 0 and particle_j.boundarydata_.boundaryparticle_==false)
        {
          double dwj=particle_j.boundarydata_.distWall_-initDist;
          double fj=1.0;
          double dmaxj=std::max(particle_j.radius_-factor_j*initDist,0.0);
          if(dwj<0.0)
            fj=0.0;
          else if(dwj<=dmaxj)
            fj=dwj/dmaxj;

          double V_j=particle_j.mass_/particle_j.density_;

          LINALG::Matrix<3,1> CFG_FF_j(particle_j.freesurfacedata_.CFG_FF_);

          sumj_Vj_CFGFFj_W_ij.Update(fj*V_j*interData_ij.w_ij_,CFG_FF_j,1.0);
          sumj_Vj_W_ij+=fj*V_j*interData_ij.w_ij_;
        }
      }//loop over j

      if(fabs(sumj_Vj_W_ij)<1.0e-10)
        dserror("Something is wrong here! To less neighbor particles, adapt factor_j!");

      // 1. step: get average CFG from surrounding particles j
      LINALG::Matrix<3,1> CFG_FF_average_i(true);
      CFG_FF_average_i.Update(1.0/sumj_Vj_W_ij,sumj_Vj_CFGFFj_W_ij,1.0);

      // 2. step: weighted sum (smooth transition) between original CFG_FF_i and extrapolated average vector CFG_FF_average_i
      LINALG::Matrix<3,1> CFG_FF_mod_i(true);
      CFG_FF_mod_i.Update(fi,CFG_FF_i,1.0);
      CFG_FF_mod_i.Update(1-fi,CFG_FF_average_i,1.0);

      double norm = CFG_FF_mod_i.Norm2();

      // 3. step: re-normalization to get a unity vector and scaling with magnitude of original CFG
      if(norm > 1.0e-12/particle_i.radius_)
        CFG_FF_mod_i.Scale(norm_CFG_FF_i/norm);
#else
      LINALG::Matrix<3,1> CFG_FF_mod_i(CFG_FF_i);
#endif

#ifdef PARTICLE_ST_BOUNDARYSCALEFAC
      //Increase the surface tension forces in the triple zone by scaling factor to account for reduced CFG magnitude due to lack of particle support
      double scale_fac=(1.0-fi)*PARTICLE_ST_BOUNDARYSCALEFAC;
      CFG_FF_mod_i.Scale(1.0+scale_fac);
#endif

      LINALG::Assemble(*CFG_FF, CFG_FF_mod_i, particle_i.lm_, particle_i.owner_);
    }// if triple point particle
    else //assemble original CFG if no triple point particle
      LINALG::Assemble(*CFG_FF, CFG_FF_i, particle_i.lm_, particle_i.owner_);

  }//loop over i

  return;
}


/*------------------------------------------------------------------------------------*
 | re-initialize density in case of free-surface flow - mesh free style  meier 09/17  |
 *------------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::MF_ReInitDensity(
    const Teuchos::RCP<Epetra_Vector> density,
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType)
{

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::MF_ReInitDensity");

  //checks
  if (density == Teuchos::null)
    dserror("density empty");

  // erase the vectors
  density->PutScalar(0.0);

  if(freeSurfaceType==INPAR::PARTICLE::InteriorReinitialization)
  {
    //loop over particles and re-initialize interior particles (FS_NONE) via density summation.
    //The density of free-surface particles has already been set via density integration and remains untouched!
    for (int lidNodeCol_i = 0; lidNodeCol_i != discret_->NodeColMap()->NumMyElements(); ++lidNodeCol_i)
    {

      ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

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

      ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

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

      ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      if(PARTICLE_REINITSHIFT==1)
      {
        if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_NONE)
          particle_i.density_=particle_i.freesurfacedata_.density_sum_;
        else
        {
          //TODO: Replace p0 with atmospheric pressue in case of inhomogeneous NBC
          double p0=0.0;

          // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
          const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);

          double density_0 = equationOfStateHandler_[matID_i]->PressureToDensity(p0,particle_i.density0_); // based on material data of fluid particle_i

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

          // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
          const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);

          double density_0 = equationOfStateHandler_[matID_i]->PressureToDensity(p0,particle_i.density0_); // based on material data of fluid particle_i

          particle_i.density_=particle_i.freesurfacedata_.density_sum_ + density_0*(1.0-particle_i.freesurfacedata_.color_field_);
        }
      }
      else if(PARTICLE_REINITSHIFT==3)
      {
        //TODO: Replace p0 with atmospheric pressue in case of inhomogeneous NBC
        double p0=0.0;

        // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
        const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);

        double density_0 = equationOfStateHandler_[matID_i]->PressureToDensity(p0,particle_i.density0_); // based on material data of fluid particle_i

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
void PARTICLE::ParticleSPHInteractionHandler::InitFreeSurfaceParticles(
    const Teuchos::RCP<Epetra_Vector> fspType)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleSPHInteractionHandler::InitFreeSurfaceParticles");

  freesurfaceparticles_.clear();

  //In the following, particles with color_field < PARTICLE_COLORLIMIT are classified as free-surface particles of type FS_DIRECT. They are summed up
  //in specific data containers freesurfaceparticles_ (containing the particles) and neighbours_fsp_ (containing their neighbours). All neighbours of
  //these particles are classified as free-surface particles of type FS_INDIRECT. These get no extra entries in freesurfaceparticles_ and neighbours_fsp_ so far!

  const int numcolelements = discret_->NodeColMap()->NumMyElements();

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeCol_i = 0; lidNodeCol_i != (unsigned int)(numcolelements); ++lidNodeCol_i)
  {
    // determine the particle_i
    ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    //Decide if norm of CFG_ is large enough to yield a reasonable approximation for the surface normal, which is not dominated by numerical noise.
    //A similar criterion is suggested by Morris et al. 2000, Eq.(20)
    if(particle_i.freesurfacedata_.CFG_FF_.Norm2()>PARTICLE_COLORGRADLIMIT/particle_i.radius_)
      particle_i.freesurfacedata_.validNormal_=true;

    if(particle_i.freesurfacedata_.color_field_ < PARTICLE_COLORLIMIT and particle_i.boundarydata_.boundaryparticle_==false)
    {
      particle_i.freesurfacedata_.freesurfaceparticletype_=FS_DIRECT;
      freesurfaceparticles_[particle_i.gid_]=&colParticles_[lidNodeCol_i];
    }
  }

  //************loop over free surface particles and search for neighbours (for free-surface particles we need superposition --> column map)********
  // bin checker
  std::set<int> examinedbins;

  // loop over the boundary particles
  for (std::map<int,ParticleSPH* >::iterator ii = freesurfaceparticles_.begin(); ii!=freesurfaceparticles_.end(); ++ii)
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
      const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

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
    ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    if(particle_i.freesurfacedata_.freesurfaceparticletype_==FS_INDIRECT)
    {
      if(particle_i.boundarydata_.boundaryparticle_==true)
        dserror("A boundary particle should not be a free-surface particle!");

      freesurfaceparticles_[particle_i.gid_]=&colParticles_[lidNodeCol_i];
    }
  }
  //***********************************************************************************************************************************************

  //************loop over free surface particles and search for neighbours (for free-surface particles we need superposition --> column map)********
  // clear bin checker
  examinedbins.clear();

  // loop over the boundary particles
  for (std::map<int,ParticleSPH* >::iterator ii = freesurfaceparticles_.begin(); ii!=freesurfaceparticles_.end(); ++ii)
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
      const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

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
      ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

      double fspType_i=particle_i.freesurfacedata_.freesurfaceparticletype_;
      LINALG::Assemble(*fspType, fspType_i, particle_i.gid_, particle_i.owner_);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | construct the specific coefficient                      sfuchs 10/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleSPHInteractionHandler::Inter_generalCoeff_ij(
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double mass_i = particle_i.mass_;

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  double density_j = 0.0;
  double mass_j = 0.0;

  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    mass_j = particle_j.mass_;
  }
  else
  {
    // In case of two-phase flow, the initial density is used for calculating the volume of a boundary particle required for generalCoeff_ij
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE");
    if(freeSurfaceType!=INPAR::PARTICLE::TwoPhase)
    {
      // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
      // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
      double pressure_j = particle_j.boundarydata_.pressureMod_;
      // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
      const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
      density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
    }
    else
    {
      density_j=particle_i.density0_;
    }
    mass_j = particle_i.mass_;
  }

  double generalCoeff_ij = 0.0;
  if(not interactionVariant2_)
    generalCoeff_ij = interData_ij.dw_ij_*mass_j;
  else
  {
    const double VSquare_ij = std::pow((mass_i/density_i), 2) + std::pow((mass_j/density_j), 2);
    generalCoeff_ij = VSquare_ij*interData_ij.dw_ij_/mass_i;
  }

  return generalCoeff_ij;
}

/*----------------------------------------------------------------------*
 | evaluate pressure                                       sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_pressure(
    LINALG::Matrix<3,1>* accn_ij,
    const double& generalCoeff_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double pressure_i = particle_i.pressure_;

  // get states of (boundary) particle j
  double density_j = 0.0;
  double pressure_j = 0.0;
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    pressure_j = particle_j.pressure_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
  }

  double fac = 0.0;
  if(not interactionVariant2_)
    fac = (pressure_i/std::pow(density_i,2) + pressure_j/std::pow(density_j,2));
  else
    fac = (density_i*pressure_j+density_j*pressure_i)/(density_i+density_j);

  accn_ij->Update(-generalCoeff_ij*fac, interData_ij.rRelVersor_ij_, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate background pressure                            sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_backgroundPressure(
    LINALG::Matrix<3,1>* trvl_accn_ij,
    const double& generalCoeff_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double pressure_i = particle_i.pressure_;
  const double radius_i = particle_i.radius_;

  // get states of (boundary) particle j
  double density_j = 0.0;
  double mass_j = 0.0;
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    mass_j = particle_j.mass_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
    mass_j = particle_i.mass_;
  }

  if(freeSurfaceType_==INPAR::PARTICLE::FreeSurface_None or freeSurfaceType_==INPAR::PARTICLE::TwoPhase)
  {
    double fac = 0.0;
    if(not interactionVariant2_)
      fac = background_pressure_*( 1.0/std::pow(density_i,2) + 1.0/std::pow(density_j,2) );
    else
      fac = background_pressure_;

    trvl_accn_ij->Update(-generalCoeff_ij*fac, interData_ij.rRelVersor_ij_, 1.0);
  }
  else
  {
    if(interactionVariant2_)
      dserror("Background pressure not available for free-surface flow in var2 so far!");

    double mod_background_pressure_i = std::min(abs(10.0*pressure_i), background_pressure_);
    double dw_ij_tilde = weightFunctionHandler_->DW(interData_ij.rRelNorm2_, radius_i/PARTICLE_P0ZHANG);
    const double p0_gradP0_Rho2_ij = mod_background_pressure_i*mass_j*dw_ij_tilde/std::pow(density_i,2);
    trvl_accn_ij->Update(-p0_gradP0_Rho2_ij, interData_ij.rRelVersor_ij_, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate xsph contribution                              sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_xsph(
    LINALG::Matrix<3,1>* trvl_accn_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double radius_i = particle_i.radius_;
  const LINALG::Matrix<3,1> vel_i(particle_i.vel_);

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  const LINALG::Matrix<3,1> vel_j( ( boundaryParticle_j ) ? particle_j.boundarydata_.velModVisc_ : particle_j.vel_ );
  // get states of (boundary) particle j
  double density_j = 0.0;
  double mass_j = 0.0;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    mass_j = particle_j.mass_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
    mass_j = particle_i.mass_;
  }

  LINALG::Matrix<3,1> vRel_ij(true);
  vRel_ij.Update(1.0, vel_i, -1.0, vel_j);

  double w_ij_damp = interData_ij.w_ij_;
  double w_ij_stiff = weightFunctionHandler_->W(interData_ij.rRelNorm2_, radius_i/PARTICLE_P0ZHANG);

  LINALG::Matrix<3,1> xsph_momentum_ij(true);
  xsph_momentum_ij.Update( -xsph_dampfac_*w_ij_damp, vRel_ij, 1.0);
  xsph_momentum_ij.Update( xsph_stiffac_*w_ij_stiff, interData_ij.rRelVersor_ij_, 1.0);

  const double inv_density_ij = 1.0/(0.5*(density_i+density_j));

  trvl_accn_ij->Update(mass_j*inv_density_ij, xsph_momentum_ij, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate additional term divA due to transport velocity sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_transportVelocity_divA(
    LINALG::Matrix<3,1>* accn_ij,
    LINALG::Matrix<3,1>* accn_ij_A,
    const double& generalCoeff_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const LINALG::Matrix<3,1> vel_i(particle_i.vel_);
  const LINALG::Matrix<3,1> velConv_i(particle_i.velConv_);

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  const LINALG::Matrix<3,1> vel_j( ( boundaryParticle_j ) ? particle_j.boundarydata_.velModVisc_ : particle_j.vel_ );
  const LINALG::Matrix<3,1> velConv_j( ( boundaryParticle_j ) ? LINALG::Matrix<3,1>(true) : particle_j.velConv_ );
  double density_j = 0.0;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
  }

  // build ternsor A_{ab}:=\rho*vel_a*(\tilde{vel}_b-vel_b) according to Adami 2013, Eq. (4).
  LINALG::Matrix<3,3> A_i(true);
  LINALG::Matrix<3,3> A_j(true);

  LINALG::Matrix<3,1> velConv_vel_i(velConv_i);
  velConv_vel_i.Update(-1.0, vel_i, 1.0);
  A_i.MultiplyNT(density_i, vel_i, velConv_vel_i, 0.0);

  // A_j is not evaluated if particle j is a boundary particle
  // Note: in the original formulation according to Adami2013, the contributions of boundary particles to the tensor A are simply set to zero
  // (according to an email by Stefan Adami, not mentioned in the paper)
  if(not boundaryParticle_j)
  {
    LINALG::Matrix<3,1> velConv_vel_j(velConv_j);
    velConv_vel_j.Update(-1.0, vel_j, 1.0);
    A_j.MultiplyNT(density_j, vel_j, velConv_vel_j, 0.0);
  }

  LINALG::Matrix<3,3> p0_A_scaled(true);
  LINALG::Matrix<3,1> p0_A_e_ij(true);

  if(not interactionVariant2_)
  {
    p0_A_scaled.Update((1.0/std::pow(density_i,2)), A_i, 1.0);
    if(not boundaryParticle_j)
      p0_A_scaled.Update((1.0/std::pow(density_j,2)), A_j, 1.0);
  }
  else
  {
    p0_A_scaled.Update(0.5, A_i, 1.0);
    if(not boundaryParticle_j)
      p0_A_scaled.Update(0.5, A_j, 1.0);
  }

  p0_A_e_ij.MultiplyNN(1.0, p0_A_scaled, interData_ij.rRelVersor_ij_, 0.0);

  if(no_veldiff_term_==false)
    accn_ij->Update(generalCoeff_ij, p0_A_e_ij, 1.0);

  accn_ij_A->Update(generalCoeff_ij, p0_A_e_ij, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate laminar viscosity                              sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_laminarViscosity(
    LINALG::Matrix<3,1>* accn_ij,
    const double& generalCoeff_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const LINALG::Matrix<3,1> vel_i(particle_i.vel_);

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  const LINALG::Matrix<3,1> vel_j( ( boundaryParticle_j ) ? particle_j.boundarydata_.velModVisc_ : particle_j.vel_ );
  double density_j = 0.0;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
  }

  // viscous interaction with boundary particles is only required in case of a no-slip boundary condition
  if(boundaryParticle_j and wallInteractionType_!=INPAR::PARTICLE::BoundarParticle_NoSlip)
    return;

  LINALG::Matrix<3,1> vRel_ij(true);
  vRel_ij.Update(1.0, vel_i, -1.0, vel_j);

  const double visc_i=particle_i.extParticleMat_->dynamicViscosity_;
  const double visc_j=particle_j.extParticleMat_->dynamicViscosity_;
  double viscosity = 0.0;
  if(visc_i>1e-12 and visc_j>1e-12)
    viscosity = 2.0*visc_i*visc_j/(visc_i+visc_j);

  if(not interactionVariant2_)
  {

    double convectionCoeff=0.0;
    double diffusionCoeff=0.0;
    const double bulkVisc_i=particle_i.extParticleMat_->bulkViscosity_;
    const double bulkVisc_j=particle_j.extParticleMat_->bulkViscosity_;
    double bulkViscosity = 0.0;
    if(bulkVisc_i>1e-12 and bulkVisc_j>1e-12)
      bulkViscosity = 2.0*bulkVisc_i*bulkVisc_j/(bulkVisc_i+bulkVisc_j);

    ViscousCoefficients(viscosity, bulkViscosity, convectionCoeff,diffusionCoeff);

    // The following two viscous terms are taken from Espanol2003, Eq(30) and re-expressed in terms of mass densities in an
    // equivalent manner. This is necessary since Espanol2003 only specifies the case m_i=m_j=m.
    const double inv_rho2rRelNorm2 = 1.0 / (density_i * density_j * interData_ij.rRelNorm2_);
    const double rRelVersorDotVrel = interData_ij.rRelVersor_ij_.Dot(vRel_ij);

    // diffusion
    const double dC_rho2rRelNorm2 = diffusionCoeff * inv_rho2rRelNorm2;
    accn_ij->Update( generalCoeff_ij*dC_rho2rRelNorm2, vRel_ij, 1.0);

    // convection
    const double cCrRelVersorDotVrel_rho2rRelNorm2 =  convectionCoeff * rRelVersorDotVrel * inv_rho2rRelNorm2;
    accn_ij->Update( generalCoeff_ij*cCrRelVersorDotVrel_rho2rRelNorm2, interData_ij.rRelVersor_ij_, 1.0);
  }
  else
  {
    accn_ij->Update( generalCoeff_ij*viscosity/interData_ij.rRelNorm2_, vRel_ij, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate artificial viscosity                           sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_artificialViscosity(
    LINALG::Matrix<3,1>* accn_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double radius_i = particle_i.radius_;
  const LINALG::Matrix<3,1> vel_i(particle_i.vel_);

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  const double radius_j = particle_j.radius_;
  const LINALG::Matrix<3,1> vel_j( ( boundaryParticle_j ) ? particle_j.boundarydata_.velModVisc_ : particle_j.vel_ );
  double density_j = 0.0;
  double mass_j = 0.0;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    mass_j = particle_j.mass_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
    mass_j = particle_i.mass_;
  }

  // viscous interaction with boundary particles is only required in case of a no-slip boundary condition
  if(boundaryParticle_j and wallInteractionType_!=INPAR::PARTICLE::BoundarParticle_NoSlip)
    return;

  LINALG::Matrix<3,1> vRel_ij(true);
  vRel_ij.Update(1.0, vel_i, -1.0, vel_j);

  // particle averaged smoothing length
  const double h_ij = 0.5*(weightFunctionHandler_->SmoothingLength(radius_i)+weightFunctionHandler_->SmoothingLength(radius_j));

  // particle averaged speed of sound
  const double c_i = particle_i.extParticleMat_->SpeedOfSoundL();
  const double c_j = particle_j.extParticleMat_->SpeedOfSoundL();
  const double c_ij = 0.5*(c_i+c_j);

  // particle averaged density
  const double density_ij = 0.5*(density_i+density_j);

  //The parameter epsilon avoids division by zero when particles get to close. The value epsilon=0.01 is a typical choice (see Adami et al. 2012).
  const double epsilon = 0.01;

  const double rRelVersorDotVrel = interData_ij.rRelVersor_ij_.Dot(vRel_ij);

  const double art_visc_fac = particle_i.extParticleMat_->artificialViscosity_*h_ij*c_ij*rRelVersorDotVrel*interData_ij.rRelNorm2_ / (density_ij*(std::pow(interData_ij.rRelNorm2_, 2)+epsilon*std::pow(h_ij, 2)));

  // the following term represents the acceleration resulting from artificial viscosity as applied in Adami et al. 2012, Eq. (11)
  accn_ij->Update(mass_j*interData_ij.dw_ij_*art_visc_fac, interData_ij.rRelVersor_ij_, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate surface tension                                sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_surfaceTension(
    LINALG::Matrix<3,1>* accn_ij,
    const ParticleSPH& particle_i,
    const ParticleSPH& particle_j,
    const InterDataPvP& interData_ij,
    const double& time
    )
{
  // get states of particle i
  const double density_i = particle_i.density_;
  const double radius_i = particle_i.radius_;
  const double mass_i = particle_i.mass_;
  // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
  const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);

  // get states of (boundary) particle j
  bool boundaryParticle_j = particle_j.boundarydata_.boundaryparticle_;
  const double radius_j = particle_j.radius_;
  double density_j = 0.0;
  double mass_j = 0.0;
  if(!boundaryParticle_j)
  {
    density_j = particle_j.density_;
    mass_j = particle_j.mass_;
  }
  else
  {
    // According to Adami et al. 2012 Eq (28), the density of boundary particle is determined based on the extrapolated pressure in Eq (27) and
    // the equation of state based on material data / mass of the interacting particle_i (for multi-phase flows, the material data might strongly differ from particle to particle).
    double pressure_j = particle_j.boundarydata_.pressureMod_;
    density_j = equationOfStateHandler_[matID_i]->PressureToDensity(pressure_j,particle_i.density0_); // based on material data of fluid particle_i
    mass_j = particle_i.mass_;
  }

  // particle averaged radius
  double radius_ij = 0.5*(radius_i+radius_j);

  // evaluate pairwise interaction potential
  double lambda = 0.0;
  double potential = 0.0;
  SurfTensionInterPot(radius_ij, interData_ij.rRelNorm2_, lambda, potential);

  if(potential != 0.0)
  {
    // particle averaged number density (n_i = \rho_i/m_i)
    double n_ij = 0.5*((density_i/mass_i)+(density_j/mass_j));

    double timefac = 1.0;
    if(time >= 0.0)
      timefac = SurfTensionTimeFac(time);

    // pair-wise interaction force scaling parameter
    double s_ij = 0.0;
    if (not boundaryParticle_j) // case of fluid-fluid interaction
    {
      if(surfaceTensionType_==INPAR::PARTICLE::ST_VDW_DIRECT)
        s_ij = 0.5*std::pow(n_ij,-2)*(surfTens_ff_/lambda);
      else
        s_ij = surfTens_ff_;
    }
    else // case of fluid-solid interaction
    {
      if(surfaceTensionType_==INPAR::PARTICLE::ST_VDW_DIRECT)
      {
        // convert static contact angle in radians
        double theta_0 = surfTens_fs_*M_PI/180.0;
        s_ij = 0.5*std::pow(n_ij,-2)*(surfTens_ff_/lambda)*(1+0.5*std::cos(theta_0));
      }
      else
        s_ij = surfTens_fs_;
    }

    // pair-wise interaction force to model surface tension following Kordilla et al. 2013 and Tartakovsky et al. 2016
    accn_ij->Update(-s_ij*potential*timefac/mass_i, interData_ij.rRelVersor_ij_, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate surface tension interaction potential          sfuchs 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::SurfTensionInterPot(
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

  // range of repulsive part
  double radius0 = surfTensPotRatio_*radius;

  if(surfaceTensionType_==INPAR::PARTICLE::ST_VDW_DIRECT)
  {
    switch (WF_DIM_)
    {
      case INPAR::PARTICLE::WF_3D :
      {
        lambda = (7.0*81.0)/(324.0*359.0) * ( -surfTensPotA_*std::pow(radius0, 2) + surfTensPotB_*std::pow(radius, 2) );
        break;
      }
      case INPAR::PARTICLE::WF_2D :
      {
        lambda = (2771.0*63.0)/(20412.0*478.0*M_PI) * ( -surfTensPotA_*std::pow(radius0, 2) + surfTensPotB_*std::pow(radius, 2) );
        break;
      }
      default :
      {
        dserror("Pairwise interaction force currently just possible for dimension 2 and 3!");
      }
    }

    if(lambda <= 0.0)
      dserror("The parameter lambda in the pairwise interaction force should be positive!");
  }

  potential = -surfTensPotA_*weightFunction->W(rRelNorm2, radius0) + surfTensPotB_*weightFunction->W(rRelNorm2, radius);

  return;
}

/*------------------------------------------------------------------------*
 | return surface tension time fac                           sfuchs 10/17 |
 *------------------------------------------------------------------------*/
double PARTICLE::ParticleSPHInteractionHandler::SurfTensionTimeFac(const double& time)
{
  double fac=1.0;
  double ramp_time=DRT::Problem::Instance()->ParticleParams().get<double>("SURFTENSION_RAMP_TIME");
  if(ramp_time>0 and time >=0)
  {
    if(time<ramp_time)
       fac=0.5*(1-cos(time*M_PI/ramp_time));
  }

  return fac;
}

/*------------------------------------------------------------------------*
 | Calculate surface forces                                   meier 10/17 |
 *------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_fspvp_Adami_1(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> curvature,
    double time)
{
  if(surfTens_ff_<=0.0)
    return;

  if(curvature!=Teuchos::null)
    curvature->PutScalar(0.0);

  // loop over the free surface particles (with superpositions -> pairs ij and ji evaluted seperately!)
  // Calculate surface curvature according to Adami et al. 2010, Eq.(20)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    //Only particles, whose color field gradient is non-zero
    if(particle_i.freesurfacedata_.CFG_FF_.Norm2()<1.0e-12/particle_i.radius_ or particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    #ifdef PARTICLE_ONLYPHASE1ACC
    if(particle_i.freesurfacedata_.phase_color_==1)
      continue;
    #endif



    double sum_nij_eij_dWij = 0.0; //numerator of Eq. (20)
    double sum_rij_dWij = 0.0; //denominator of Eq. (20)

    //denominator according to Morris et al. 2000, Eq. (24). Attention: Here we need a special color field and cannot reuse the color field already calculated earlier,
    //since here only the particle with a sufficiently high color field gradient are considered!!!
    double sum_Wij = 0.0;
    //auto contribution
    sum_Wij+=weightFunctionHandler_->W0(particle_i.radius_) * particle_i.mass_ /particle_i.density_;

    //unity normal vector
    LINALG::Matrix<3,1> n_i(particle_i.freesurfacedata_.CFG_F2S_);
    double norm_n_i=n_i.Norm2();
    n_i.Scale(1.0/norm_n_i);

    // loop over the interaction particle list
    boost::unordered_map<int, InterDataPvP>::const_iterator jj;
    for (jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // determine particle j
      const int& lidNodeCol_j = jj->first;
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij = jj->second;

      //Only particles, whose color field gradient is non-zero
      if(particle_j.freesurfacedata_.CFG_FF_.Norm2()<1.0e-12/particle_i.radius_  or particle_j.freesurfacedata_.CFG_F2S_.Norm2()<1.0e-10*particle_i.radius_ or particle_j.boundarydata_.boundaryparticle_==true)
        continue;

      if (interData_ij.dw_ij_ != 0)
      {
        //unity normal vector
        LINALG::Matrix<3,1> n_j(particle_j.freesurfacedata_.CFG_F2S_);
        double norm_n_j=n_j.Norm2();
        n_j.Scale(1.0/norm_n_j);

        double V_j=particle_j.mass_/particle_j.density_;
        sum_nij_eij_dWij+=(n_i.Dot(interData_ij.rRelVersor_ij_)-n_j.Dot(interData_ij.rRelVersor_ij_))*V_j*interData_ij.dw_ij_;
        sum_rij_dWij+=interData_ij.rRelNorm2_*V_j*interData_ij.dw_ij_;
        sum_Wij+=V_j*interData_ij.w_ij_;
      }
    }

    double kappa_i=0.0;
    if(abs(sum_Wij)>1.0e-10)
      kappa_i=-sum_nij_eij_dWij/sum_Wij;

    //TODO: Alternative normalization according to Adami et al. 2010 (uncomment if required)
//    double dim=(double)(WF_DIM_)+1.0;
//    if(abs(sum_rij_dWij)>1.0e-10)
//      kappa_i=dim*sum_nij_eij_dWij/sum_rij_dWij;

    LINALG::Matrix<3,1> acc_i(particle_i.freesurfacedata_.CFG_FF_);

    if(surfTens_ff_>0.0)
    {
      double timefac = 1.0;
      if(time >= 0.0)
        timefac = SurfTensionTimeFac(time);

      acc_i.Scale(-timefac*kappa_i*surfTens_ff_/particle_i.density_);

      LINALG::Assemble(*curvature, kappa_i, particle_i.gid_, particle_i.owner_);
      LINALG::Assemble(*accn, acc_i, particle_i.lm_, particle_i.owner_);
    }

  }

  return;
}

/*------------------------------------------------------------------------*
 | Calculate surface forces                                   meier 10/17 |
 *------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_fspvp_Adami_2(
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> curvature,
    double time)
{
  if(surfTens_ff_<=0.0)
    return;

  if(curvature!=Teuchos::null)
    curvature->PutScalar(0.0);

  // loop over the free surface particles (with superpositions -> pairs ij and ji evaluted seperately!)
  // Calculate surface curvature according to Adami et al. 2010, Eq.(20)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    //Only particles, whose color field gradient is large enough (in order to avoid the noise of small contributions). --> See Morris et al. 2000, Eq.(20)
    if(particle_i.freesurfacedata_.validNormal_==false or particle_i.boundarydata_.boundaryparticle_==true)
      continue;

    double sum_nij_eij_dWij = 0.0; //numerator of Eq. (20)
    double sum_rij_dWij = 0.0; //denominator of Eq. (20)

    //denominator according to Morris et al. 2000, Eq. (24). Attention: Here we need a special color field and cannot reuse the color field already calculated earlier,
    //since here only the particle with a sufficiently high color field gradient are considered!!!
    double sum_Wij = 0.0;
    //auto contribution
    sum_Wij+=weightFunctionHandler_->W0(particle_i.radius_) * particle_i.mass_ /particle_i.density_;

    //unity normal vector
    LINALG::Matrix<3,1> n_i(particle_i.freesurfacedata_.CFG_FF_);
    double norm_n_i=n_i.Norm2();
    n_i.Scale(1.0/norm_n_i);

    // loop over the interaction particle list
    boost::unordered_map<int, InterDataPvP>::const_iterator jj;
    for (jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // determine particle j
      const int& lidNodeCol_j = jj->first;
      const ParticleSPH& particle_j = colParticles_[lidNodeCol_j];
      const InterDataPvP& interData_ij = jj->second;

    //Only particles, whose color field gradient is large enough (in order to avoid the noise of small contributions). --> See Morris et al. 2000, Eq.(20)
    if(particle_j.freesurfacedata_.validNormal_==false or particle_j.boundarydata_.boundaryparticle_==true)
      continue;

      if (interData_ij.dw_ij_ != 0)
      {
        //unity normal vector
        LINALG::Matrix<3,1> n_j(particle_j.freesurfacedata_.CFG_FF_);
        double norm_n_j=n_j.Norm2();
        n_j.Scale(1.0/norm_n_j);

        double V_j=particle_j.mass_/particle_j.density_;
        sum_nij_eij_dWij+=(n_i.Dot(interData_ij.rRelVersor_ij_)-n_j.Dot(interData_ij.rRelVersor_ij_))*V_j*interData_ij.dw_ij_;
        sum_rij_dWij+=interData_ij.rRelNorm2_*V_j*interData_ij.dw_ij_;
        sum_Wij+=V_j*interData_ij.w_ij_;
      }
    }

    double kappa_i=0.0;
    //TODO
    if(abs(sum_Wij)>1.0e-10)
      kappa_i=-sum_nij_eij_dWij/sum_Wij;

    //TODO: Alternative normalization according to Adami et al. 2010 (uncomment if required)
//    double dim=(double)(WF_DIM_)+1.0;
//    if(abs(sum_rij_dWij)>1.0e-10)
//      kappa_i=dim*sum_nij_eij_dWij/sum_rij_dWij;

    LINALG::Matrix<3,1> acc_i(particle_i.freesurfacedata_.CFG_FF_);
    double expfac=1.0;

#ifdef PARTICLE_SFEXP
    expfac=1.0/std::pow(abs(particle_i.freesurfacedata_.color_field_),PARTICLE_SFEXP);
#endif

    if(surfTens_ff_>0.0)
    {
      double timefac = 1.0;
      if(time >= 0.0)
        timefac = SurfTensionTimeFac(time);

      acc_i.Scale(-expfac*timefac*kappa_i*surfTens_ff_/particle_i.density_);

      LINALG::Assemble(*curvature, kappa_i, particle_i.gid_, particle_i.owner_);
      LINALG::Assemble(*accn, acc_i, particle_i.lm_, particle_i.owner_);
    }

  }

  return;
}

/*-------------------------------------------------------------------*
 Alternative model to enforce the static contact angle   meier 10/17 |
 *-------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Inter_fspvp_Adami_ContactAngleVar2(
    const Teuchos::RCP<Epetra_Vector> accn,
    double time)
{
  if(surfTens_ff_<=0.0)
    return;

  // loop over the free surface particles (with superpositions -> pairs ij and ji evaluted seperately!)
  // Calculate surface curvature according to Adami et al. 2010, Eq.(20)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    //Only particles, whose color field gradient is non-zero
    if(particle_i.freesurfacedata_.CFG_FF_.Norm2()<1.0e-12/particle_i.radius_ or particle_i.freesurfacedata_.CFG_F2S_.Norm2()<1.0e-10*particle_i.radius_
        or particle_i.boundarydata_.boundaryparticle_==true or particle_i.boundarydata_.haveboundaryneighbor_==false)
      continue;

    #ifdef PARTICLE_ONLYPHASE1ACC
    if(particle_i.freesurfacedata_.phase_color_==1)
      continue;
    #endif

    LINALG::Matrix<3,1> normal1(particle_i.freesurfacedata_.CFG_FF_);
    double norm1=normal1.Norm2();

    if(norm1 > 1.0e-10*particle_i.radius_)
      normal1.Scale(1.0/norm1);
    else
      normal1.Scale(0.0);

    LINALG::Matrix<3,1> normal_wi(particle_i.freesurfacedata_.CFG_F1S_);
    double normal_wi_norm=normal_wi.Norm2();
    if(normal_wi_norm>1.0e-10*particle_i.radius_)
      normal_wi.Scale(1.0/normal_wi_norm);
    else
      normal_wi.Scale(0.0);

    LINALG::Matrix<3,1> normal_ti(true);
    normal_ti.Update(1.0,normal1,1.0);
    normal_ti.Update(-normal_wi.Dot(normal1),normal_wi,1.0);
    double norm_ti=normal_ti.Norm2();
    if(norm_ti > 1.0e-10*particle_i.radius_)
      normal_ti.Scale(1.0/norm_ti);
    else
      normal_ti.Scale(0.0);

    double cosalpha=-normal_wi.Dot(normal1);
    double cosalpha0=cos(surfTens_fs_*M_PI/180.0);

    //The factor 2.0 is required since the boundary dirac delta functioin is only evaluated in the fluid and not in the solid phase (integral gives 0.5 and not 1)!
    double force=surfTens_ff_*(cosalpha-cosalpha0)*particle_i.freesurfacedata_.CFG_FF_.Norm2()*2.0*normal_wi_norm;

    LINALG::Matrix<3,1> acc_i(normal_ti);

    double timefac = 1.0;

    if(time >= 0.0)
      timefac = SurfTensionTimeFac(time);

    acc_i.Scale(timefac*force/particle_i.density_);

    LINALG::Assemble(*accn, acc_i, particle_i.lm_, particle_i.owner_);
  }

  return;
}

/*------------------------------------------------------------------------------*
 | Determine coefficients required for viscous forces              meier 10/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::ViscousCoefficients(const double viscosity, const double bulkViscosity, double& convectionCoeff, double& diffusionCoeff)
{

#ifndef PARTICLE_ONLYLAPLACETERM
//Attention: The trace of a tensor as required to derive equation (28) from equation (25) in Espanol2003 depends
//on the spatial dimension --> SPH approximations of the Laplace operator typically differ for different dimensions!
switch (WF_DIM_)
{
  case INPAR::PARTICLE::WF_3D :
  {
    convectionCoeff = 5.0*(bulkViscosity+viscosity / 3.0);
    break;
  }
  case INPAR::PARTICLE::WF_2D :
  {
    convectionCoeff = 4.0*(bulkViscosity+viscosity / 3.0);
    break;
  }
  case INPAR::PARTICLE::WF_1D :
  {
    convectionCoeff = 3.0*(bulkViscosity+viscosity / 3.0);
    break;
  }
  default :
  {
    dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }
}

diffusionCoeff = 5.0 * viscosity / 3.0 - bulkViscosity;
#else
diffusionCoeff = 2.0 * viscosity;
convectionCoeff = 0.0;
if(bulkViscosity>0)
  dserror("Bulk Viscosity not considered here!");
#endif

// checks
if (diffusionCoeff<0)
{
  dserror("The diffusion coefficient is negative! The following equation should hold: 5*dynamicViscosity >= 3*bulkViscosity");
}

if (convectionCoeff<0)
{
  dserror("The convection coefficient is negative! Are you sure that the dynamic viscosity and the bulk modulus are positive?");
}

  return;
}

/*------------------------------------------------------------------------------*
 | Compute pressure vector (see e.g. Antoci2007-E4)                meier 10/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::Density2Pressure(
  const Teuchos::RCP<const Epetra_Vector> density,
  Teuchos::RCP<Epetra_Vector> &pressure)
{
  // checks
  if (density == Teuchos::null)
  {
    pressure = Teuchos::null;
    return;
  }

  // build the pressure vector
  if (pressure == Teuchos::null)
    pressure = Teuchos::rcp(new Epetra_Vector(density->Map(), true));

  // compute inertia for every particle
  for (int lidNodeRow_i = 0; lidNodeRow_i < density->MyLength(); ++lidNodeRow_i)
  {
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    (*pressure)[lidNodeRow_i] = equationOfStateHandler_[matID_i]->DensityToPressure((*density)[lidNodeRow_i],particle_i.density0_);
  }
}

/*------------------------------------------------------------------------------*
 | Determine internal energy                                       meier 10/17  |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleSPHInteractionHandler::DetermineIntEnergy(double &intenergy)
{
  intenergy=0.0;
  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    // determine particle i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleSPH& particle_i = colParticles_[lidNodeCol_i];

    // Boundary particles are not considered here (there is no potential existent for the pressure / density values extrapolated to the boundary particles)
    // On the long term, a more elaborate scheme might be desirable that at least consideres the work contributions of boundary particles!
    if(particle_i.boundarydata_.boundaryparticle_)
      continue;

    // determine phase particle_i belongs to and material (boundary particles have phase color -1 and get the same material as phase 0)
    const int matID_i = std::max(particle_i.freesurfacedata_.phase_color_,0);
    intenergy += equationOfStateHandler_[matID_i]->DensityToEnergy(particle_i.density_,particle_i.mass_,particle_i.density0_);
  }
  return;
}
