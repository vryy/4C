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
#include "../drt_adapter/ad_str_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
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
  const double restDensity) :
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  myrank_(discret->Comm().MyPID()),
  weightFunctionHandler_(Teuchos::null),
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE")),
  restDensity_(restDensity),
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
  }
}

// other checks
if (wallInteractionType_ != INPAR::PARTICLE::Mirror)
{
  if (wallMeshFreeData_.density_ < 0)
    dserror("the value of WALL_FAKE_DENSITY is unacceptable");
  if (wallMeshFreeData_.mass_ < 0)
    dserror("the value of WALL_FAKE_MASS is unacceptable");
}

//Attention: The trace of a tensor as required to derive equation (28) from equation (25) in Espanol2003 depends
//on the spatial dimension --> SPH approximations of the Laplace operator typically differ for differnt dimensions!
switch (PARTICLE_DIM)
{
  case 3 :
  {
    diffusionCoeff_ = 5.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
    break;
  }
  case 2 :
  {
    diffusionCoeff_ = 8.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
    break;
  }
  case 1 :
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
  weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_CubicBspline);
  break;
}
case INPAR::PARTICLE::SqrtHyperbola :
{
  weightFunctionHandler_ = Teuchos::rcp(new WeightFunction_SqrtHyperbola);
  break;
}
case INPAR::PARTICLE::HyperbolaNoRsz :
{
  weightFunctionHandler_ = Teuchos::rcp(new WeightFunction_HyperbolaNoRsz);
  break;
}
}

}


/*----------------------------------------------------------------------*
 | set up internal variables for future computations       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Init(
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> densityn,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> pressure,
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<const Epetra_Vector> densityapproxn,
    const int step)
{
  // check
  if (colParticles_.size() != 0)
  {
    dserror("you did not call Clear before Init, colParticles_ is not empty");
  }
  // set up the local data storage and fill it with the state vectors
  InitColParticles(disn, veln, radiusn, mass, densityn, specEnthalpyn, pressure, temperature, densityapproxn);

  // set neighbours
  AddNewNeighbours(step);

  // security block
  trg_updatedColorFieldGradient_ = false;

  //Printneighbours_p();
}


/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::InitColParticles(
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> densityn,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> pressure,
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<const Epetra_Vector> densityapproxn)
{
  // row to col vectors
    // dof-based vectors
  Teuchos::RCP<Epetra_Vector> disnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
  Teuchos::RCP<Epetra_Vector> velnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
    // node-based vectors
  Teuchos::RCP<Epetra_Vector> radiusnCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> massCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> densitynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> specEnthalpynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> pressureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> temperatureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
  Teuchos::RCP<Epetra_Vector> densityapproxCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
    // exports
  LINALG::Export(*disn ,*disnCol);
  LINALG::Export(*veln ,*velnCol);
  LINALG::Export(*radiusn ,*radiusnCol);
  LINALG::Export(*mass ,*massCol);
  LINALG::Export(*densityn ,*densitynCol);
  LINALG::Export(*specEnthalpyn ,*specEnthalpynCol);
  LINALG::Export(*pressure ,*pressureCol);
  LINALG::Export(*temperature ,*temperatureCol);
  LINALG::Export(*densityapproxn ,*densityapproxCol);

  const int numcolelements = discret_->NodeColMap()->NumMyElements();

  colParticles_.resize(numcolelements);
  for (int lidNodeCol=0; lidNodeCol<numcolelements; ++lidNodeCol)
  {
    DRT::Node *particle = discret_->lColNode(lidNodeCol);
    std::vector<int> lm;
    lm.reserve(3);
    discret_->Dof(particle, lm);

    // dof-based vectors
    LINALG::Matrix<3,1> dis, vel;
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disnCol, dis, lm);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velnCol, vel, lm);
    // node-based vectors
    const double radius = (*radiusnCol)[lidNodeCol];
    const double mass = (*massCol)[lidNodeCol];
    const double density = (*densitynCol)[lidNodeCol];
    const double specEnthalpy = (*specEnthalpynCol)[lidNodeCol];
    const double pressure = (*pressureCol)[lidNodeCol];
    const double temperature = (*temperatureCol)[lidNodeCol];
    const double densityapprox = (*densityapproxCol)[lidNodeCol];
    // set up the particles
    colParticles_[lidNodeCol] = ParticleMF(particle->Id(), particle->Owner(), lm,
        dis, vel, radius, mass, density, specEnthalpy, pressure, temperature, densityapprox);
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
  }
}
}


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

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn(Teuchos::null);
  Teuchos::RCP<const Epetra_Vector> wallveln(Teuchos::null);
  if(walldiscret != Teuchos::null)
  {
    walldisn = walldiscret->GetState("walldisnp");
    wallveln = walldiscret->GetState("wallvelnp");
  }

  // bin checker
  std::set<int> examinedbins;

  // loop over the particles (no superpositions)
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
    particle_algorithm_->GetNeighbouringItems(binId, neighboursLinf_p, neighboursLinf_w, neighboursLinf_hs);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMF& particle_i = colParticles_[lidNodeCol_i];

      AddNewNeighbours_p(particle_i, neighboursLinf_p, step);

      AddNewNeighbours_w(particle_i, neighboursLinf_w, walldiscret, walldisn, wallveln, step);

      AddNewNeighbours_hs(particle_i, neighboursLinf_hs);
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

  // loop over the neighbours_ particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboursLinf_p.begin(); jj!=neighboursLinf_p.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

    // evaluate contact only once
    if(particle_i.gid_ < particle_j.gid_)
    {
      // create the data that we have to insert
      LINALG::Matrix<3,1> rRel_ij;
      rRel_ij.Update(1.0, particle_i.dis_, -1.0, particle_j.dis_); // inward vector
      const double rRelNorm2 = rRel_ij.Norm2();

      const double weightDerivative_ij = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_i.radius_);
      const double weight_ij = weightFunctionHandler_->Weight(rRelNorm2, particle_i.radius_);
      double weightDerivative_ji = 0;
      double weight_ji = 0;
      if (particle_j.owner_ == myrank_)
      {
        if (particle_i.radius_ == particle_j.radius_)
        {
          weightDerivative_ji = weightDerivative_ij;
          weight_ji = weight_ij;
        }
        else
        {
          weightDerivative_ji = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_j.radius_);
          weight_ji = weightFunctionHandler_->Weight(rRelNorm2, particle_j.radius_);
        }
      }

      // push_back
      if (weightDerivative_ij != 0 || weightDerivative_ji != 0)
      {
        LINALG::Matrix<3,1> rRelVersor_ij(rRel_ij);
        rRelVersor_ij.Scale(1/rRelNorm2);

        const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
        (neighbours_p_[lidNodeRow_i])[lidNodeCol_j] = InterDataPvP(
            rRelVersor_ij,
            weightDerivative_ij,
            weightDerivative_ji,
            weight_ij,
            weight_ji,
            rRelNorm2,
            step);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - walls                              katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::AddNewNeighbours_w(
    const ParticleMF& particle_i,
    const boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
    const Teuchos::RCP<DRT::Discretization>& walldiscret,
    const Teuchos::RCP<const Epetra_Vector>& walldisn,
    const Teuchos::RCP<const Epetra_Vector>& wallveln,
    const int step)
{
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
    const double weightDerivative = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_i.radius_);

    // insert in the set if appropriate
    if (weightDerivative != 0)
    {
      LINALG::Matrix<3,1> rRelVersor(rRel);
      rRelVersor.Scale(1/rRelNorm2);

      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid_);
      InterDataPvW& interData_ij = (neighbours_w_[lidNodeRow_i])[wallcontact.elemId_];
      interData_ij = InterDataPvW(
        rRelVersor,
        weightDerivative,
        rRelNorm2,
        wallcontact.point_,
        step);
    }
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
 | update weights and distances in all the neighbours - particles     katta 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights_p(const int step)
{
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

        const double weightDerivative_ij = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_i.radius_);
        const double weight_ij = weightFunctionHandler_->Weight(rRelNorm2, particle_i.radius_);
        double weightDerivative_ji = 0;
        double weight_ji = 0;
        if (particle_j.owner_ == myrank_)
        {
          if (particle_i.radius_ == particle_j.radius_)
          {
            weightDerivative_ji = weightDerivative_ij;
            weight_ji = weight_ij;
          }
          else
          {
            weightDerivative_ji = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_j.radius_);
            weight_ji = weightFunctionHandler_->Weight(rRelNorm2, particle_j.radius_);
          }
        }

        // insert the new weights, do not touch the step
        interData_ij = InterDataPvP(
          rRelVersor_ij,
          weightDerivative_ij,
          weightDerivative_ji,
          weight_ij,
          weight_ji,
          rRelNorm2,
          interData_ij.step_);
      }
    }
  }
}


/*---------------------------------------------------------------------------------*
 | update weights and distances in all the neighbours - particles     katta 01/17  |
 *---------------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::UpdateWeights_w(const int step)
{
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
      const int& lidNodeCol_j = jj->first;
      InterDataPvW& interData_ij =(neighbours_w_[lidNodeRow_i])[lidNodeCol_j];
      if (step != interData_ij.step_)
      {
        // --- extract and compute general data --- //

        // create the data that we have to insert
        LINALG::Matrix<3,1> rRel;
        rRel.Update(1.0, particle_i.dis_, -1.0, interData_ij.contactPoint_); // inward vector
        const double rRelNorm2 = rRel.Norm2();
        LINALG::Matrix<3,1> rRelVersor(rRel);
        rRelVersor.Scale(1/rRelNorm2);
        const double weightDerivative = weightFunctionHandler_->WeightDerivative(rRelNorm2, particle_i.radius_);

        // insert the new weights, do not touch the step
        interData_ij = InterDataPvW(
            rRelVersor,
            weightDerivative,
            rRelNorm2,
            interData_ij.contactPoint_,
            interData_ij.step_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | print interactions                                      katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Print_neighbours_p()
{
  std::cout << "\n\n particle - particle interactions\n\n";
  //bool trg_interactions = false;
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours_p_.size(); ++lidNodeRow_i)
  {
    std::cout << discret_->NodeRowMap()->GID(lidNodeRow_i) << " |";

    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      std::cout << " " << colParticles_[jj->first].gid_;
      //trg_interactions = true;
    }
    std::cout << std::endl;
  }
}


/*----------------------------------------------------------------------*
 | evaluate density - pvp                                  katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_densityDot(
    const Teuchos::RCP<Epetra_Vector> densityDotn)
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
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // --- extract and compute general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);
      const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);
      const double density = rRelVersorDotVrel;

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
        // assemble and write
        double densityDotn_ij = generalCoeff_ij * density;
        LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid_, particle_i.owner_);

      }

      // write on particle j if appropriate specializing the quantities
      if ( interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        double densityDotn_ji = generalCoeff_ji * density;
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
    const bool trg_pressure)
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

    // determine some useful quantities
    const double rhoSquare_i = std::pow(particle_i.density_,2);
    const double p_Rho2_i = particle_i.pressure_/rhoSquare_i;

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      LINALG::Matrix<3,1> vRel_ij;
      vRel_ij.Update(1.0, particle_i.vel_, -1.0, particle_j.vel_);
      const double rho2 = particle_i.density_ * particle_j.density_;
      const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

      LINALG::Matrix<3,1> momentum_ij;
        // pressure
      if (trg_pressure)
      {
        const double rhoSquare_j = std::pow(particle_j.density_,2);
        const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure_/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
        momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);
      }
        // diffusion
      const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
      momentum_ij.Update(dC_rho2rRelNorm2, vRel_ij, 1.0);
      // convection
      const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
      momentum_ij.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor_ij, 1.0);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        LINALG::Assemble(*accn, accn_ij, particle_i.lm_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        accn_ji.Scale(-1.0); // actio = - reactio
        LINALG::Assemble(*accn, accn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate specEnthalpyDot - pvp                          katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_specEnthalpyDot(
    const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
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

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // --- extract general data --- //

      const double rRelNorm2 = interData_ij.rRelNorm2_;
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      const double rho2 = particle_i.density_ * particle_j.density_;
      // assemble and write
      const double deltaT_ij = particle_i.temperature_ - particle_j.temperature_;
      const double divT_rho_ij = 2 * particle_algorithm_->ExtParticleMat()->thermalConductivity_ *
          deltaT_ij / (rho2 * rRelNorm2);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;

        // --- specific enthalpy --- //

        double specEnthalpyDotn_ij = generalCoeff_ij * divT_rho_ij;
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        double specEnthalpyDotn_ji = generalCoeff_ji * divT_rho_ij;
        specEnthalpyDotn_ji *= -1.0; // actio = - reactio
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ji, particle_j.gid_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate color field gradient - pvp                     katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_colorFieldGradient(
    const Teuchos::RCP<Epetra_Vector> colorFieldGradientn)
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

      // --- extract general data --- //
      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);

      LINALG::Matrix<3,1> colorFieldGradientGeneral_ij;
      colorFieldGradientGeneral_ij.Update(particle_i.radius_ / particle_j.density_,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> colorFieldGradientn_ij;
        colorFieldGradientn_ij.Update(generalCoeff_ij, colorFieldGradientGeneral_ij);
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ij, particle_i.lm_, particle_i.owner_);

      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> colorFieldGradientn_ji;
        colorFieldGradientn_ji.Update(generalCoeff_ji, colorFieldGradientGeneral_ij);
        colorFieldGradientn_ji.Scale(-1.0); // actio = - reactio
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_surfaceTensionCFG(
  const Teuchos::RCP<Epetra_Vector> accn)
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

    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // --- extract general data --- //

      const double densityCorrectiveTerm = 2 * restDensity_ / (particle_i.density_ + particle_j.density_);

      // rescaling and assembpling
      LINALG::Matrix<3,1> acc_ij;
      acc_ij.Update(1.0, particle_i.colorFieldGradient_, -1.0, particle_j.colorFieldGradient_);
      acc_ij.Scale(- densityCorrectiveTerm * surfaceVoidTension_);
      LINALG::Assemble(*accn, acc_ij, particle_i.lm_, particle_i.owner_);

      // rescaling and assembpling
      LINALG::Matrix<3,1> acc_ji;
      acc_ji.Update(-1.0,acc_ij); // actio = - reaction
      LINALG::Assemble(*accn, acc_ji, particle_j.lm_, particle_j.owner_);
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate divergence free pressures - pvp                katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvp_divFreePressureAcc(
    const Teuchos::RCP<Epetra_Vector> divFreePressureAcc,
    const double dt)
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

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      const double kv_rho_j = particle_j.densityDot_ * particle_j.alpha_ / (particle_j.density_ * dt);

      const double kv_rho = - (kv_rho_i + kv_rho_j);
      LINALG::Matrix<3,1> momentum_ij;
      momentum_ij.Update(kv_rho,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        LINALG::Assemble(*divFreePressureAcc, accn_ij, particle_i.lm_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        accn_ji.Scale(-1.0); // actio = - reactio
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
    const double dt)
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
    const double k_rho_i = (particle_i.density_ - restDensity_) * particle_i.alpha_ / (particle_i.density_ * dt * dt);

    // loop over the interaction particle list
    for (boost::unordered_map<int, InterDataPvP>::const_iterator jj = neighbours_p_[lidNodeRow_i].begin(); jj != neighbours_p_[lidNodeRow_i].end(); ++jj)
    {
      // extract data for faster access
      const int& lidNodeCol_j = jj->first;
      const InterDataPvP& interData_ij = jj->second;
      const ParticleMF& particle_j = colParticles_[lidNodeCol_j];

      // --- extract general data --- //

      LINALG::Matrix<3,1> rRelVersor_ij(interData_ij.rRelVersor_ij_);
      const double k_rho_j = (particle_j.density_ - restDensity_) * particle_j.alpha_ / (particle_j.density_ * dt * dt);

      const double kv_rho = - (k_rho_i + k_rho_j);
      LINALG::Matrix<3,1> momentum_ij;
      momentum_ij.Update(kv_rho,rRelVersor_ij);

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
        const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ij;
        accn_ij.Update(generalCoeff_ij, momentum_ij);
        LINALG::Assemble(*constDensityPressureAcc, accn_ij, particle_i.lm_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ji coeff
        const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
        // assemble and write
        LINALG::Matrix<3,1> accn_ji;
        accn_ji.Update(generalCoeff_ji, momentum_ij);
        accn_ji.Scale(-1.0); // actio = - reactio
        LINALG::Assemble(*constDensityPressureAcc, accn_ji, particle_j.lm_, particle_j.owner_);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate density - pvw                                  katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_densityDot(
    const Teuchos::RCP<Epetra_Vector> densityDotn)
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

      // construct the specific coeff
      const double generalCoeff = interData_ij.weightDerivative_ * wallMeshFreeData_.mass_;
      // assemble and write
      double densityDotn_ij = generalCoeff * density;
      LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid_, particle_i.owner_);
    }
  }
}


/*----------------------------------------------------------------------*
 | evaluate acceleration - pvw                             katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvw_acc(
    const Teuchos::RCP<Epetra_Vector> accn)
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

    // determine some useful quantities
    const double rhoSquare_i = std::pow(particle_i.density_,2);
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
         wallMeshFreeData_.density_ = particle_i.density_;
         wallMeshFreeData_.mass_ = particle_i.mass_;
         // mirrored velocity
         vRel.Scale(2);
       }

       const double rho2 = particle_i.density_ * wallMeshFreeData_.density_;
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

       // construct the specific coeff
       const double generalCoeff = interData_ij.weightDerivative_ * wallMeshFreeData_.mass_;
       // assemble and write
       LINALG::Matrix<3,1> accn_ij;
       accn_ij.Update(generalCoeff, momentum);
       LINALG::Assemble(*accn, accn_ij, particle_i.lm_, particle_i.owner_);
    }
  }
}



/*--------------------------------------------------------------------------*
 | evaluate specEnthalpyDot - pvhs                             katta 10/16  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Inter_pvhs_specEnthalpyDot(
    const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
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


  double specEnthalpyDot_i = 0.0;
  // loop over the interaction particle list
  for (boost::unordered_map<int, Teuchos::RCP<HeatSource> >::const_iterator jj = neighbours_hs_[lidNodeRow_i].begin(); jj != neighbours_hs_[lidNodeRow_i].end(); ++jj)
  {
    const Teuchos::RCP<HeatSource> hs = jj->second;
    specEnthalpyDot_i += (hs->QDot_)/particle_i.density_;
  }
  LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDot_i, particle_i.gid_, particle_i.owner_);
}
}


/*--------------------------------------------------------------------------*
 | compute density - mesh free style                           katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_density(const Teuchos::RCP<Epetra_Vector> densityn)
{
  //checks
  if (densityn == Teuchos::null)
  {
    dserror("densityn is empty");
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

      // write on particle i if appropriate specializing the quantities
      if (interData_ij.weight_ij_ != 0)
      {
        // assemble and write
        double density_ij = interData_ij.weight_ij_ * particle_j.mass_;
        LINALG::Assemble(*densityn, density_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weight_ji_ != 0)
      {
        // assemble and write
        double density_ji = interData_ij.weight_ji_ * particle_i.mass_;
        LINALG::Assemble(*densityn, density_ji, particle_j.gid_, particle_j.owner_);
      }
    }
  }
}


/*--------------------------------------------------------------------------*
 | compute density - mesh free style                           katta 01/17  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::MF_alpha()
{
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
      if (interData_ij.weightDerivative_ij_ != 0)
      {
        // construct the specific ij coeff
         const double generalCoeff_ij = interData_ij.weightDerivative_ij_ * particle_j.mass_;
         // assemble and write dof
         LINALG::Matrix<3,1> alphaDof_ij;
         alphaDof_ij.Update(generalCoeff_ij, rRelVersor_ij);
         LINALG::Assemble(*alphaDof, alphaDof_ij, particle_i.lm_, particle_i.owner_);
        // assemble and write node
        double alpha_ij = alphaDof_ij.Dot(alphaDof_ij);
        LINALG::Assemble(*alpha, alpha_ij, particle_i.gid_, particle_i.owner_);
      }

      // write on particle j if appropriate specializing the quantities
      if (interData_ij.weightDerivative_ji_ != 0)
      {
        // construct the specific ij coeff
         const double generalCoeff_ji = interData_ij.weightDerivative_ji_ * particle_i.mass_;
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
  SetStateVector(alpha, Alpha);
}

