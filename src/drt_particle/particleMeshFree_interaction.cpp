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
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE")),
  restDensity_(restDensity),
  trg_secondRound_(false)
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
    wallMeshFreeData_.density = extParticleMat->initDensity_;
    wallMeshFreeData_.mass = extParticleMat->initDensity_ * PARTICLE::Utils::Radius2Volume(extParticleMat->initRadius_);
    // the pressure is linked to the deltaDensity_ with the initial density. In case of this wall it is always 0
    wallMeshFreeData_.pressure = 0;
    /*
    if (extParticleMat->initTemperature_ < extParticleMat->transitionTemperature_)
    {
      wallMeshFreeData_.pressure = PARTICLE::Utils::Density2Pressure(extParticleMat->SpeedOfSoundS(),extParticleMat->initDensity_);
    }
    else if (extParticleMat->initTemperature_ > extParticleMat->transitionTemperature_)
    {
      wallMeshFreeData_.pressure = PARTICLE::Utils::Density2Pressure(extParticleMat->SpeedOfSoundL(),extParticleMat->initDensity_);
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
    wallMeshFreeData_.density = -1;
    wallMeshFreeData_.mass = -1;
    wallMeshFreeData_.pressure = -1;
    break;
  }
case INPAR::PARTICLE::Custom :
  {
    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    wallMeshFreeData_.density = particleparams.get<double>("WALL_FAKE_DENSITY");
    wallMeshFreeData_.mass = particleparams.get<double>("WALL_FAKE_MASS");
    wallMeshFreeData_.pressure = particleparams.get<double>("WALL_FAKE_PRESSURE");
  }
}

// other checks
if (wallInteractionType_ != INPAR::PARTICLE::Mirror)
{
  if (wallMeshFreeData_.density < 0)
    dserror("the value of WALL_FAKE_DENSITY is unacceptable");
  if (wallMeshFreeData_.mass < 0)
    dserror("the value of WALL_FAKE_MASS is unacceptable");
}

diffusionCoeff_ = 5.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkViscosity_;
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
  weightDerivative_ = &PARTICLE::WeightFunction_CubicBspline::WeightDerivative;
  break;
}
case INPAR::PARTICLE::SqrtHyperbola :
{
  weightDerivative_ = &PARTICLE::WeightFunction_SqrtHyperbola::WeightDerivative;
  break;
}
case INPAR::PARTICLE::HyperbolaNoRsz :
{
  weightDerivative_ = &PARTICLE::WeightFunction_HyperbolaNoRsz::WeightDerivative;
  break;
}
}

// do we need the second round? here we decide
if (surfaceVoidTension_ != 0)
{
  trg_secondRound_ = true;
}


}


/*----------------------------------------------------------------------*
 | set up internal variables for future computations       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Init(
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> densityn,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<const Epetra_Vector> pressure)
{
  if (particleMeshFreeData_.size() == 0 && neighbouringParticles_.size() ==0)
  {
    // set state vectors
    SetStateVectors(disn,veln,radiusn,densityn,specEnthalpyn,mass,temperature,pressure);

    // set neighbours
    SetNeighbours();
  }
  else
  {
    dserror("The temporary content of the particle interactions is not empty");
  }

  //PrintNeighbouringParticles();
}


/*----------------------------------------------------------------------*
 | clear data, keep memory                                 katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::Clear()
{
  // erase particleMeshFreeData_. keep the memory
  particleMeshFreeData_.clear();
  // erase neighbours keep memory
  neighbouringParticles_.clear();
  neighbouringWalls_.clear();
  neighbouringHeatSources_.clear();
}

/*----------------------------------------------------------------------*
 | set colVectors in the local data structs                katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVectors(
    Teuchos::RCP<const Epetra_Vector> disn,
    Teuchos::RCP<const Epetra_Vector> veln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<const Epetra_Vector> densityn,
    Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
    Teuchos::RCP<const Epetra_Vector> mass,
    Teuchos::RCP<const Epetra_Vector> temperature,
    Teuchos::RCP<const Epetra_Vector> pressure)
{
// checks
if (disn == Teuchos::null     ||
    veln == Teuchos::null     ||
    radiusn == Teuchos::null  ||
    densityn == Teuchos::null ||
    mass == Teuchos::null     ||
    pressure == Teuchos::null ||
    temperature == Teuchos::null)
{
  dserror("one or more state vectors are empty");
}

/// miraculous transformation into column vectors... ///

// dof based vectors
Teuchos::RCP<Epetra_Vector> disnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*disn,*disnCol);
Teuchos::RCP<Epetra_Vector> velnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*veln,*velnCol);
// node based vectors
Teuchos::RCP<Epetra_Vector> radiusnCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*radiusn,*radiusnCol);
Teuchos::RCP<Epetra_Vector> densitynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*densityn,*densitynCol);
Teuchos::RCP<Epetra_Vector> specEnthalpynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*specEnthalpyn,*specEnthalpynCol);
Teuchos::RCP<Epetra_Vector> massCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*mass,*massCol);
Teuchos::RCP<Epetra_Vector> temperatureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*temperature,*temperatureCol);
Teuchos::RCP<Epetra_Vector> pressureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*pressure,*pressureCol);

// fill particleData_
const int numcolparticles = discret_->NodeColMap()->NumMyElements();
particleMeshFreeData_.resize(numcolparticles);

for (int i=0; i<numcolparticles; ++i)
{
  // particle for which data will be collected
  DRT::Node *particle = discret_->lColNode(i);

  std::vector<int> lm;
  lm.reserve(3);

  // extract global dof ids and fill into lm_i
  discret_->Dof(particle, lm);

  ParticleMeshFreeData& data = particleMeshFreeData_[particle->LID()];

  //position, velocity and angular velocity of particle
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disnCol,data.dis,lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velnCol,data.vel,lm);

  const int lid = particle->LID();
  // node-based state vectors
  data.gid = particle->Id();
  data.radius = (*radiusnCol)[lid];
  data.density = (*densitynCol)[lid];
  data.specEnthalpy = (*specEnthalpynCol)[lid];
  data.mass = (*massCol)[lid];
  data.pressure = (*pressureCol)[lid];
  data.temperature = (*temperatureCol)[lid];


  // lm vector and owner
  data.lm.swap(lm);
  data.owner = particle->Owner();

#ifdef DEBUG
if(particle->NumElement() != 1)
  dserror("More than one element for this particle");
#endif
}

}

/*----------------------------------------------------------------------*
 | set all the neighbours                                  katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetNeighbours()
{
  // resize the vectors
  const int numRowParticles = discret_->NodeRowMap()->NumMyElements();
  neighbouringParticles_.resize(numRowParticles);
  neighbouringWalls_.resize(numRowParticles);
  neighbouringHeatSources_.resize(numRowParticles);

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

    // list of walls that border on the CurrentBin
    std::set<DRT::Element*> neighboring_walls;

    std::list<DRT::Node*> neighboring_particles;

    // list of heat sources that border on the CurrentBin
    const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources = Teuchos::rcp(new std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less>);

    // first neighbouring round
    particle_algorithm_->GetNeighbouringItems(binId, neighboring_particles, neighboring_walls, neighboring_heatSources);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing the the important addresses
      const int lidNodeCol_i = currentBinParticles[i]->LID();
      const ParticleMeshFreeData& particle_i = particleMeshFreeData_[lidNodeCol_i];


      SetNeighbours_Particles(particle_i, neighboring_particles);

      SetNeighbours_Walls(particle_i, neighboring_walls, walldiscret, walldisn, wallveln);

      SetNeighbours_HeatSources(particle_i, *neighboring_heatSources);
    }
  }
}


/*----------------------------------------------------------------------*
 | set the neighbours - particles                          katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetNeighbours_Particles(
    const ParticleMeshFreeData& particle_i,
    const std::list<DRT::Node*>& neighboring_particles)
{
  // self-neighbours not allowed
  // insert the interaction only if meaningful

  // loop over the neighbouring particles
  for(std::list<DRT::Node*>::const_iterator jj=neighboring_particles.begin(); jj!=neighboring_particles.end(); ++jj)
  {
    const int lidNodeCol_j = (*jj)->LID();
    const ParticleMeshFreeData &particle_j = particleMeshFreeData_[lidNodeCol_j];

    // evaluate contact only once
    if(particle_i.gid < particle_j.gid)
    {

      // create the data that we have to push_back
      Interaction_Particles data;

      data.lidCol_j = lidNodeCol_j;
      data.rRelVersor_ij.Update(1.0, particle_i.dis, -1.0, particle_j.dis); // inward vector
      data.rRelNorm2 = data.rRelVersor_ij.Norm2();
      data.rRelVersor_ij.Scale(1/data.rRelNorm2);

      // compute the proper weight function gradient
      data.weightDerivative_ij = weightDerivative_(data.rRelNorm2, particle_i.radius);
      if (particle_j.owner == myrank_)
      {
        if (particle_i.radius == particle_j.radius)
        {
          data.weightDerivative_ji = data.weightDerivative_ij;
        }
        else
        {
          data.weightDerivative_ji = weightDerivative_(data.rRelNorm2, particle_j.radius);
        }
      }
      else
      {
        data.weightDerivative_ji = 0;
      }

      // push_back
      if (data.weightDerivative_ij != 0 || data.weightDerivative_ji != 0)
      {
        const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid);
        neighbouringParticles_[lidNodeRow_i].push_back(data);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | set the neighbours - walls                              katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetNeighbours_Walls(
    const ParticleMeshFreeData& particle_i,
    const std::set<DRT::Element*>& neighboring_walls,
    const Teuchos::RCP<DRT::Discretization>& walldiscret,
    const Teuchos::RCP<const Epetra_Vector>& walldisn,
    const Teuchos::RCP<const Epetra_Vector>& wallveln)
{
  // evaluate contact with walls first
  std::vector<WallInteractionPoint> surfaces;
  std::vector<WallInteractionPoint> lines;
  std::vector<WallInteractionPoint> nodes;

  std::set<int> unusedIds;

  // check whether there is contact between particle i and neighboring walls
  for(std::set<DRT::Element*>::const_iterator w=neighboring_walls.begin(); w!=neighboring_walls.end();  ++w)
  {
    DRT::Element* neighboringwallele = (*w);
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
    GEO::ObjectType objecttype = GEO::nearest3DObjectOnElement(neighboringwallele,nodeCoord,particle_i.dis,nearestPoint);
    //-----------------------------------------------------------------------------------------

    static LINALG::Matrix<3,1> r_i_wall;
    r_i_wall.Update(1.0, nearestPoint, -1.0, particle_i.dis);
    const double distance_i_wall = r_i_wall.Norm2();
    const double penetration = distance_i_wall-particle_i.radius;

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
        distance_vector.Update(1.0, nearestPoint, -1.0, (*pointer)[i].point);
        const double distance = distance_vector.Norm2();
        const double adaptedtol = GEO::TOL7 * particle_i.radius;

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
        WallInteractionPoint currentContact = { neighboringwallele->Id(), nearestPoint, penetration, nodeCoord, lm_wall, lmowner };
        (*pointer).push_back(currentContact);
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
    const double rminusg = particle_i.radius -std::abs(surfaces[s].penetration);
    const double radius_surface = sqrt(particle_i.radius * particle_i.radius - rminusg * rminusg);

    for(size_t l=0; l<lines.size(); ++l)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, surfaces[s].point, -1.0, lines[l].point);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_surface)
        unusedIds.insert(lines[l].elemId);
    }
    for(size_t p=0; p<nodes.size(); ++p)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, surfaces[s].point, -1.0, nodes[p].point);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_surface)
        unusedIds.insert(nodes[p].elemId);
    }
  }
  // find entries of nodes which are within the penetration volume of the current particle
  // hierarchical: lines next
  for(size_t l=0; l<lines.size(); ++l)
  {
    // radius = sqrt(r_i^2 - (r_i-|g|)^2)
    const double rminusg = particle_i.radius - std::abs(lines[l].penetration);
    const double radius_line = sqrt(particle_i.radius * particle_i.radius - rminusg*rminusg);

    for(size_t p=0; p<nodes.size(); ++p)
    {
      static LINALG::Matrix<3,1> distance_vector;
      distance_vector.Update(1.0, lines[l].point, -1.0, nodes[p].point);
      const double distance = distance_vector.Norm2();
      if(distance <= radius_line)
        unusedIds.insert(nodes[p].elemId);
    }
  }

  // write entries of lines and nodes to surfaces if contact has to be evaluated
  for(size_t l=0; l<lines.size(); ++l)
    if( !unusedIds.count(lines[l].elemId) )
      surfaces.push_back(lines[l]);
  for(size_t p=0; p<nodes.size(); ++p)
    if( !unusedIds.count(nodes[p].elemId) )
      surfaces.push_back(nodes[p]);

  // evaluate contact between particle_i and entries of surfaces
  const int lidNodeCol_i = discret_->NodeColMap()->LID(particle_i.gid);
  std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(discret_->lColNode(lidNodeCol_i))->Get_history_wall();
  if(history_wall.size() > 3)
    dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    WallInteractionPoint wallcontact = surfaces[s];

    // create the data that we have to push_back
    Interaction_Walls data;

    data.elemId = wallcontact.elemId;
    data.rRelVersor.Update(1.0, particle_i.dis, -1.0, wallcontact.point); // inward vector
    data.rRelNorm2 = data.rRelVersor.Norm2();
    data.rRelVersor.Scale(1/data.rRelNorm2);

    // compute the proper weight function gradient
    data.weightDerivative = weightDerivative_(data.rRelNorm2, particle_i.radius);

    // push_back
    if (data.weightDerivative != 0)
    {
      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid);
      neighbouringWalls_[lidNodeRow_i].push_back(data);
    }
  }
}



/*----------------------------------------------------------------------*
 | set the neighbours - heat sources                       katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::SetNeighbours_HeatSources(
    const ParticleMeshFreeData& particle_i,
    const std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less>& neighboring_heatSources)
{
  std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less>::const_iterator hs;
  for(hs = neighboring_heatSources.begin(); hs != neighboring_heatSources.end();  ++hs)
  {
    if ((*hs)->minVerZone_[0]<=particle_i.dis(0) &&
        (*hs)->minVerZone_[1]<=particle_i.dis(1) &&
        (*hs)->minVerZone_[2]<=particle_i.dis(2) &&
        (*hs)->maxVerZone_[0]>=particle_i.dis(0) &&
        (*hs)->maxVerZone_[1]>=particle_i.dis(1) &&
        (*hs)->maxVerZone_[2]>=particle_i.dis(2))
    {
      const int lidNodeRow_i = discret_->NodeRowMap()->LID(particle_i.gid);
      neighbouringHeatSources_[lidNodeRow_i].push_back(*hs);
    }
  }
}


/*----------------------------------------------------------------------*
 | print interactions                                      katta 12/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::PrintNeighbouringParticles()
{
  std::cout << "particle - particle interactions\n\n";
  //bool trg_interactions = false;
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbouringParticles_.size(); ++lidNodeRow_i)
  {
    std::cout << discret_->NodeRowMap()->GID(lidNodeRow_i) << " |";

    for (std::list<Interaction_Particles>::const_iterator jj = neighbouringParticles_[lidNodeRow_i].begin(); jj != neighbouringParticles_[lidNodeRow_i].end(); ++jj)
    {
      std::cout << " " << particleMeshFreeData_[jj->lidCol_j].gid;
      //trg_interactions = true;
    }
    std::cout << std::endl;
  }
}


/*----------------------------------------------------------------------*
 | set colVectors in the local data structs                katta 10/16  |
 *----------------------------------------------------------------------*/
// Between E2 and E3 in http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf
// The input must be a row vector
void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVectors_SecondRound(
    Teuchos::RCP<const Epetra_Vector> colorFieldGradientn)
{
// checks
if (colorFieldGradientn == Teuchos::null)
{
  dserror("one or more state vectors are empty");
}

/// miraculous transformation into column vectors... ///

// dof based vectors
Teuchos::RCP<Epetra_Vector> colorFieldGradientnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*colorFieldGradientn,*colorFieldGradientnCol);

// fill particleData_
const int numcolparticles = discret_->NodeColMap()->NumMyElements();

  for (int i=0; i<numcolparticles; ++i)
  {
    // particle for which data will be collected
    DRT::Node *particle = discret_->lColNode(i);

    ParticleMeshFreeData& data = particleMeshFreeData_[particle->LID()];

    //colorFieldGradientCol of particle
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*colorFieldGradientnCol,data.colorFieldGradient,data.lm);
  }
}

/*----------------------------------------------------------------------*
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::EvaluateParticleMeshFreeInteractions(
  const Teuchos::RCP<Epetra_Vector> accn,
  const Teuchos::RCP<Epetra_Vector> densityDotn,
  const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
{
  //checks
  if (accn == Teuchos::null || densityDotn == Teuchos::null || specEnthalpyDotn == Teuchos::null)
  {
    dserror("one or more input vectors are empty");
  }

  Teuchos::RCP<Epetra_Vector> colorFieldGradientn = Teuchos::null;
  if (trg_secondRound_)
  {
    colorFieldGradientn = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }

  // loop over the particles (no superpositions)
  for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbouringParticles_.size(); ++lidNodeRow_i)
  {
    // determine the particle_i
    const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
    const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
    const ParticleMeshFreeData& particle_i = particleMeshFreeData_[lidNodeCol_i];

    CalcNeighboringHeatSourcesContact(particle_i, neighbouringHeatSources_[lidNodeRow_i], specEnthalpyDotn);

    CalcNeighboringParticleMeshFreeInteraction(particle_i, neighbouringParticles_[lidNodeRow_i], accn, densityDotn, specEnthalpyDotn, colorFieldGradientn);

    CalcNeighboringWallMeshFreeInteraction(particle_i, neighbouringWalls_[lidNodeRow_i], accn, densityDotn);
  }


  // second round
  if (trg_secondRound_)
  {
    // reorder the interesting vectors
    SetStateVectors_SecondRound(colorFieldGradientn);

    // loop over the particles (no superpositions)
    for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbouringParticles_.size(); ++lidNodeRow_i)
    {
      // determine the particle_i
      const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);
      const int lidNodeCol_i = discret_->NodeColMap()->LID(gid_i);
      const ParticleMeshFreeData& particle_i = particleMeshFreeData_[lidNodeCol_i];

      CalcNeighboringParticleMeshFreeInteraction_SecondRound(particle_i, neighbouringParticles_[lidNodeRow_i], accn);

    }
  }
}


/*----------------------------------------------------------------------*
 | calc neighbouring particleMeshFree interactions         katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringParticleMeshFreeInteraction(
  const ParticleMeshFreeData& particle_i,
  const std::list<Interaction_Particles> neighbouringParticles,
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<Epetra_Vector> densityDotn,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
  Teuchos::RCP<Epetra_Vector> colorFieldGradientn)
{
  // check, is the list empty?
  if (neighbouringParticles.size() == 0)
  {
    return;
  }

  // determine some useful quantities
  const double rhoSquare_i = std::pow(particle_i.density,2);
  const double p_Rho2_i = particle_i.pressure/rhoSquare_i;

  // loop over the interaction particle list
  for (std::list<Interaction_Particles>::const_iterator jj = neighbouringParticles.begin(); jj != neighbouringParticles.end(); ++jj)
  {
    const Interaction_Particles& interactionData = *jj;

    const ParticleMeshFreeData& particle_j = particleMeshFreeData_[interactionData.lidCol_j];

    // --- extract general data --- //

    const double rRelNorm2 = interactionData.rRelNorm2;
    LINALG::Matrix<3,1> rRelVersor_ij(interactionData.rRelVersor_ij);
    LINALG::Matrix<3,1> vRel_ij;
    vRel_ij.Update(1.0, particle_i.vel, -1.0, particle_j.vel);
    const double rho2 = particle_i.density * particle_j.density;
    const double rRelVersorDotVrel = rRelVersor_ij.Dot(vRel_ij);

    // --- density --- //
    const double density = rRelVersorDotVrel;

    // --- momentum --- //
    LINALG::Matrix<3,1> momentum_ij;
      // pressure
    const double rhoSquare_j = std::pow(particle_j.density,2);
    const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
    momentum_ij.Update(- gradP_Rho2, rRelVersor_ij, 1.0);
      // diffusion
    const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
    momentum_ij.Update(dC_rho2rRelNorm2, vRel_ij, 1.0);
    // convection
    const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
    momentum_ij.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor_ij, 1.0);
    // correction parameter for adhesion-cohesion mechanics
    // const double densityCorrectiveTerm = 2 * restDensity_ / (particle_i.density + particle_j.density);
    // compute the cohesion term
    //const double cohesionWeight = PARTICLE::SurfaceTensionInteractions::Cohesion(rRelNorm2, particle_i.radius);
    //momentum_i.Update( - densityCorrectiveTerm * surfaceVoidTension_ * cohesionWeight,rRelVersor,1.0);

    // --- specific enthalpy --- //
    const double deltaT_ij = particle_i.temperature - particle_j.temperature;
    const double divT_rho_ij = 2 * particle_algorithm_->ExtParticleMat()->thermalConductivity_ *
        deltaT_ij / (rho2 * rRelNorm2);

    // --- color field --- //
    LINALG::Matrix<3,1> colorFieldGradientGeneral_ij;
    if (trg_secondRound_)
    {
      colorFieldGradientGeneral_ij.Update(particle_i.radius / particle_j.density,rRelVersor_ij);
    }

    // write on particle i if appropriate specializing the quantities
    if (interactionData.weightDerivative_ij != 0)
    {
      // construct the specific ij coeff
      const double generalCoeff_ij = interactionData.weightDerivative_ij * particle_j.mass;

      // --- density --- //

      double densityDotn_ij = generalCoeff_ij * density;
      LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid, particle_i.owner);

      // --- acceleration --- //

      LINALG::Matrix<3,1> accn_ij;
      accn_ij.Update(generalCoeff_ij, momentum_ij);
      LINALG::Assemble(*accn, accn_ij, particle_i.lm, particle_i.owner);

      // --- specific enthalpy --- //

      double specEnthalpyDotn_ij = generalCoeff_ij * divT_rho_ij;
      LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ij, particle_i.gid, particle_i.owner);

      // --- color field --- //
      if (trg_secondRound_)
      {
        // mass scaling and assembling
        LINALG::Matrix<3,1> colorFieldGradientn_ij;
        colorFieldGradientn_ij.Update(generalCoeff_ij, colorFieldGradientGeneral_ij);
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ij, particle_i.lm, particle_i.owner);
      }
    }

    // write on particle j if appropriate specializing the quantities
    if (interactionData.weightDerivative_ji != 0)
    {
      // construct the specific ji coeff
      const double generalCoeff_ji = interactionData.weightDerivative_ji * particle_i.mass;

      // --- density --- //

      double densityDotn_ji = generalCoeff_ji * density;
      LINALG::Assemble(*densityDotn, densityDotn_ji, particle_j.gid, particle_j.owner);

      // --- acceleration --- //

      LINALG::Matrix<3,1> accn_ji;
      accn_ji.Update(generalCoeff_ji, momentum_ij);
      accn_ji.Scale(-1.0); // actio = - reactio
      LINALG::Assemble(*accn, accn_ji, particle_j.lm, particle_j.owner);

      // --- specific enthalpy --- //

      double specEnthalpyDotn_ji = generalCoeff_ji * divT_rho_ij;
      specEnthalpyDotn_ji *= -1.0; // actio = - reactio
      LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_ji, particle_j.gid, particle_j.owner);

      // --- color field --- //
      if (trg_secondRound_)
      {
        // mass scaling and assembling
        LINALG::Matrix<3,1> colorFieldGradientn_ji;
        colorFieldGradientn_ji.Update(generalCoeff_ji, colorFieldGradientGeneral_ij);
        colorFieldGradientn_ji.Scale(-1.0); // actio = - reactio
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_ji, particle_j.lm, particle_j.owner);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | compute interaction with neighboring walls              katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringWallMeshFreeInteraction(
    const ParticleMeshFreeData& particle_i,
    const std::list<Interaction_Walls>& neighbouringWalls,
    const Teuchos::RCP<Epetra_Vector> accn,
    const Teuchos::RCP<Epetra_Vector> densityDotn)
{
  // check, is the list empty?
  if (neighbouringWalls.size() == 0)
  {
    return;
  }

  const double rhoSquare_i = std::pow(particle_i.density,2);
  const double p_Rho2_i = particle_i.pressure/rhoSquare_i;

  for (std::list<Interaction_Walls>::const_iterator jj = neighbouringWalls.begin(); jj != neighbouringWalls.end(); ++jj)
  {
    const Interaction_Walls& interactionData = *jj;

    // --- extract general data --- //

    const double rRelNorm2 = interactionData.rRelNorm2;
    LINALG::Matrix<3,1> rRelVersor(interactionData.rRelVersor);
    LINALG::Matrix<3,1> vRel(particle_i.vel);

    if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
    {
      wallMeshFreeData_.pressure = particle_i.pressure;
      wallMeshFreeData_.density = particle_i.density;
      wallMeshFreeData_.mass = particle_i.mass;
      // mirrored velocity
      vRel.Scale(2);
    }

    const double rho2 = particle_i.density * wallMeshFreeData_.density;
    const double rRelVersorDotVrel = rRelVersor.Dot(vRel);

    // --- density --- //
    const double density = rRelVersorDotVrel;

    // --- momentum --- //
    LINALG::Matrix<3,1> momentum;
      // pressure
    const double rhoSquare_j = std::pow(wallMeshFreeData_.density,2);
    const double gradP_Rho2 = (p_Rho2_i + wallMeshFreeData_.pressure/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
    momentum.Update(- gradP_Rho2, rRelVersor, 1.0);
      // diffusion
    const double dC_rho2rRelNorm2 = diffusionCoeff_ / (rho2 * rRelNorm2);
    momentum.Update(dC_rho2rRelNorm2, vRel, 1.0);
    // convection
    const double cCrRelVersorDotVrel_rho2rRelNorm2 = convectionCoeff_ * rRelVersorDotVrel / (rho2 * rRelNorm2);
    momentum.Update(cCrRelVersorDotVrel_rho2rRelNorm2, rRelVersor, 1.0);

    // construct the specific coeff
    const double generalCoeff = interactionData.weightDerivative * wallMeshFreeData_.mass;

    // --- density --- //

    double densityDotn_ij = generalCoeff * density;
    LINALG::Assemble(*densityDotn, densityDotn_ij, particle_i.gid, particle_i.owner);

    // --- acceleration --- //

    LINALG::Matrix<3,1> accn_ij;
    accn_ij.Update(generalCoeff, momentum);
    LINALG::Assemble(*accn, accn_ij, particle_i.lm, particle_i.owner);

  }
}

/*--------------------------------------------------------------------------*
 | calculate interaction with neighboring heat sources         katta 10/16  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringHeatSourcesContact(
    const ParticleMeshFreeData& particle_i,
  const std::list<Teuchos::RCP<HeatSource> > neighboringHeatSources,
  const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn)
{
  double specEnthalpyDot_i = 0.0;

  std::list<Teuchos::RCP<HeatSource> >::const_iterator hs;
  for(hs = neighboringHeatSources.begin(); hs != neighboringHeatSources.end();  ++hs)
  {
      specEnthalpyDot_i += ((*hs)->QDot_)/particle_i.density;
  }

  LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDot_i, particle_i.gid, particle_i.owner);
}

/*----------------------------------------------------------------------*
 | calc neighbouring particleMeshFree interactions         katta 10/16  |
 *----------------------------------------------------------------------*/
// no check of the distance, this is a second round without weightDerivative!
void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringParticleMeshFreeInteraction_SecondRound(
  const ParticleMeshFreeData &particle_i,
  const std::list<Interaction_Particles> neighbouringParticles,
  const Teuchos::RCP<Epetra_Vector> accn)
{
  // check, is the list empty?
  if (neighbouringParticles.size() == 0)
  {
    return;
  }

  // loop over the interaction particle list
  for (std::list<Interaction_Particles>::const_iterator jj = neighbouringParticles.begin(); jj != neighbouringParticles.end(); ++jj)
  {
    const Interaction_Particles& interactionData = *jj;

    const ParticleMeshFreeData& particle_j = particleMeshFreeData_[interactionData.lidCol_j];

    // --- extract general data --- //

    const double densityCorrectiveTerm = 2 * restDensity_ / (particle_i.density + particle_j.density);

    // rescaling and assembpling
    LINALG::Matrix<3,1> acc_ij;
    acc_ij.Update(1.0, particle_i.colorFieldGradient, -1.0, particle_j.colorFieldGradient);
    acc_ij.Scale(- densityCorrectiveTerm * surfaceVoidTension_);

    // --- acceleration --- //

    LINALG::Assemble(*accn, acc_ij, particle_i.lm, particle_i.owner);

    // rescaling and assembpling
    LINALG::Matrix<3,1> acc_ji;
    acc_ji.Update(-1.0,acc_ij); // actio = - reaction

    // --- acceleration --- //

    LINALG::Assemble(*accn, acc_ji, particle_j.lm, particle_j.owner);
  }
}
