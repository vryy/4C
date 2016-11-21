/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.cpp

\brief Particle-MeshFree interaction handling

\level 3

\maintainer Alessandro Cattabiani
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particleMeshFree_interaction.H"
#include "particleMeshFree_weightFunction.H"
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
  const Teuchos::ParameterList& particledynparams) :
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  myrank_(discret->Comm().MyPID()),
  weightFunctionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION")),
  wallInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE"))
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

diffusionCoeff_ = 5.0 * extParticleMat->dynamicViscosity_ / 3.0 - extParticleMat->bulkModulus_;
convectionCoeff_ = 5.0 * (extParticleMat->dynamicViscosity_ / 3.0  + extParticleMat->bulkModulus_);

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
    Teuchos::RCP<const Epetra_Vector> pressure)
{
// checks
if (disn == Teuchos::null     ||
    veln == Teuchos::null     ||
    radiusn == Teuchos::null  ||
    densityn == Teuchos::null ||
    mass == Teuchos::null     ||
    pressure == Teuchos::null)
  dserror("one or more state vectors are empty");

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
Teuchos::RCP<Epetra_Vector> pressureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*pressure,*pressureCol);
// create the temperature row vector
Teuchos::RCP<const Epetra_Vector> temperature = PARTICLE::Utils::SpecEnthalpy2Temperature(specEnthalpyn,particle_algorithm_->ExtParticleMat());
Teuchos::RCP<Epetra_Vector> temperatureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*temperature,*temperatureCol);

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
 | set colVectors in the local data structs                katta 10/16  |
 *----------------------------------------------------------------------*/
/// --- second round -> unactivated ----------------------------------- ///
/*
void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVectors_SecondRound(
    Teuchos::RCP<const Epetra_Vector> gradT_over_rho)
{
// checks
if (gradT_over_rho == Teuchos::null)
  dserror("one or more state vectors are empty");

/// miraculous transformation into column vectors... ///

// dof based vectors
Teuchos::RCP<Epetra_Vector> gradT_over_rhoCol = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*gradT_over_rho,*gradT_over_rhoCol);

// fill particleData_
const int numcolparticles = discret_->NodeColMap()->NumMyElements();

  for (int i=0; i<numcolparticles; ++i)
  {
    // particle for which data will be collected
    DRT::Node *particle = discret_->lColNode(i);

    ParticleMeshFreeData& data = particleMeshFreeData_[particle->LID()];

    //gradT_over_rhoCol of particle
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*gradT_over_rhoCol,data.gradT_over_rho,data.lm);

  }
}
*/
///--------------------------------------------------------------------///

/*----------------------------------------------------------------------*
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::EvaluateParticleMeshFreeInteractions(
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<Epetra_Vector> densityDotn,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
{
  //checks
  if (accn == Teuchos::null || densityDotn == Teuchos::null)
    dserror("one or more input vectors are empty");

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn(Teuchos::null);
  Teuchos::RCP<const Epetra_Vector> wallveln(Teuchos::null);
  if(walldiscret != Teuchos::null)
  {
    walldisn = walldiscret->GetState("walldisnp");
    wallveln = walldiscret->GetState("wallvelnp");
  }

  // store bins, which have already been examined
  std::set<int> examinedbins;

/// --- second round -> unactivated ----------------------------------- ///
//  // vector to collect the first part of the specEnthalpy computations
//  Teuchos::RCP<Epetra_Vector> gradT_over_rho = LINALG::CreateVector(*discret_->DofRowMap(),true);
//
//  // list of the list of all particles in the neighborhood of currparticle and bins for multiple iterations
//  std::list<DRT::Element*> binList;
//  std::list<std::list<DRT::Node*> > neighboring_particlesList;
///--------------------------------------------------------------------///

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

    particle_algorithm_->GetNeighbouringItems(binId, neighboring_particles, neighboring_walls, neighboring_heatSources);

/// --- second round -> unactivated ----------------------------------- ///
//    // store for future use
//    binList.push_back(currentBin);
//    neighboring_particlesList.push_back(neighboring_particles);
///--------------------------------------------------------------------///

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing
      DRT::Node* particle_i = currentBinParticles[i];
      const int lidNodeCol_i = particle_i->LID();

      // compute neighbouring heat sources
      CalcNeighboringHeatSourcesContact(lidNodeCol_i, neighboring_heatSources, specEnthalpyDotn);

      // compute contact with neighboring walls
      CalcNeighboringWallMeshFreeInteraction(lidNodeCol_i, neighboring_walls, walldiscret, walldisn, wallveln, accn);

      // compute interactions with neighboring particles
      CalcNeighboringParticleMeshFreeInteraction(lidNodeCol_i, neighboring_particles, accn, densityDotn, specEnthalpyDotn);


    }
  }


/// --- second round -> unactivated ----------------------------------- ///
//  // reorder the interesting vectors for the second round
//  SetStateVectors_SecondRound(gradT_over_rho);
//
//  std::list<std::list<DRT::Node*> >::const_iterator currNeighboring_particles = neighboring_particlesList.begin();
//  for (std::list<DRT::Element*>::const_iterator currBin = binList.begin(); currBin != binList.end(); ++currBin)
//  {
//    // extract the pointer to the particles
//    DRT::Node** currentBinParticles = (*currBin)->Nodes();
//
//    // loop over all particles in CurrentBin
//    for(int i=0; i<(*currBin)->NumNode(); ++i)
//    {
//      // determine the particle we are analizing
//      DRT::Node* particle_i = currentBinParticles[i];
//      const int lidNodeCol_i = particle_i->LID();
//
//      // compute interactions with neighboring particles - second round
//      CalcNeighboringParticleMeshFreeInteraction_SecondRound(lidNodeCol_i, (*currNeighboring_particles), specEnthalpyDotn);
//
//    }
//
//    ++currNeighboring_particles;
//  }
///--------------------------------------------------------------------///

  // erase temporary storage for interaction data. keep the memory
  particleMeshFreeData_.clear();
}


/*----------------------------------------------------------------------*
 | calc neighbouring particleMeshFree interactions         katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringParticleMeshFreeInteraction(
  const int lidNodeCol_i,
  const std::list<DRT::Node*> neighboring_particles,
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<Epetra_Vector> densityDotn,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
{
  // ids i
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];
  const double rhoSquare_i = std::pow(particle_i.density,2);
  const double p_Rho2_i = particle_i.pressure/rhoSquare_i;
  //const double T_over_rhoSquare_i = particle_i.temperature/rhoSquare_i;

  // self-interaction
  // densityDot -> the weightFunction gradient is null (or ill-posed in case of a strange weightFunction)
  // acc -> the weightFunction gradient is null (or ill-posed in case of a strange weightFunction)

  // loop over the neighbouring particles
  // avoid to recompute forces (self-interaction are not allowed in this part of the code)
  for(std::list<DRT::Node*>::const_iterator jj=neighboring_particles.begin(); jj!=neighboring_particles.end(); ++jj)
  {
    const ParticleMeshFreeData &particle_j = particleMeshFreeData_[(*jj)->LID()];

    // evaluate contact only once in case we own particle j and the radii match. Otherwise compute everything,
    // the assemble method does not write in case the particle is a ghost.
    // another check on the radii is performed (later in the code) to handle the case where we want to write, we can write but accelerations are not actio<->reactio
    if(particle_i.gid < particle_j.gid || particle_i.radius != particle_j.radius || particle_j.owner != myrank_)
    {

      // compute distance and relative velocities
      LINALG::Matrix<3,1> rRel, vRel;
      double WFDerivative = 0;
      rRel.Update(1.0, particle_i.dis, -1.0, particle_j.dis);
      const double rRelNorm2 = rRel.Norm2();
      vRel.Update(1.0, particle_i.vel, -1.0, particle_j.vel);
      const double rRelDotVrel = rRel.Dot(vRel);

      // compute the proper weight function gradient
      switch (weightFunctionType_)
      {
      case INPAR::PARTICLE::CubicBspline :
      {
        WFDerivative = PARTICLE::WeightFunction_CubicBspline::DerivativeWeight(rRelNorm2, particle_i.radius);
        break;
      }
      }

      // useful quantities
      const double rhoSquare_j = std::pow(particle_j.density,2);
      // p_i/\rho_i^2 + p_j/\rho_j^2
      const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure/rhoSquare_j);
      // WFDerivative/(\rho_i \rho_j)
      const double WFDerivative_Rho2 = WFDerivative / (particle_i.density * particle_j.density);
      // T_i/\rho_i^2 - T_j/\rho_j^2
      //const double gradT_over_rho_coeff = (particle_j.temperature/rhoSquare_j - T_over_rhoSquare_i);

      // compute WFGradDotVrel
      const double WFGradDotVrel = WFDerivative * rRelDotVrel;

      // compute divT_rho_i
      const double divT_rho_i = 2 * particle_algorithm_->ExtParticleMat()->thermalConductivity_ * WFDerivative_Rho2 * (particle_i.temperature - particle_j.temperature);

      LINALG::Matrix<3,1> momentum_i;
      // compute the pressure term
      momentum_i.Update(- gradP_Rho2 * WFDerivative,rRel);
      // compute the diffusion term
      momentum_i.Update(diffusionCoeff_ * WFDerivative_Rho2,vRel,1.0);
      // compute the convection term
      momentum_i.Update(convectionCoeff_ * WFDerivative_Rho2 * rRelDotVrel / (rRelNorm2 * rRelNorm2),rRel,1.0);

      // mass scalings
      double densityDotn_i = particle_j.mass * WFGradDotVrel;
      double specEnthalpyDotn_i = particle_j.mass * divT_rho_i;
      LINALG::Matrix<3,1> accn_i;
      accn_i.Update(particle_j.mass,momentum_i);

      // compute and add the correct densityDot_i
      LINALG::Assemble(*densityDotn, densityDotn_i, particle_i.gid, particle_i.owner);
      // compute and add the correct specEnthalpyDot_i
      LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_i, particle_i.gid, particle_i.owner);
      // compute and add the correct accn_i
      LINALG::Assemble(*accn, accn_i, particle_i.lm, particle_i.owner);
      // compute and add the correct gradT_over_rho_i
      //LINALG::Assemble(*gradT_over_rho, gradT_over_rho_i, particle_i.lm, particle_i.owner);


      // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
      // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
      if(particle_i.radius == particle_j.radius)
      {
        double densityDotn_j = particle_i.mass * WFGradDotVrel; // actio = - reactio (but twice! nothing changes)
        double specEnthalpyDotn_j = - particle_i.mass * divT_rho_i; // actio = - reactio -> divT_rho_j = - divT_rho_i
        LINALG::Matrix<3,1> accn_j;//, gradT_over_rho_j(WFGrad);
        accn_j.Update(- particle_i.mass,momentum_i); // actio = - reactio -> momentum_j = - momentum_i
        //gradT_over_rho_j.Scale(particle_i.mass * gradT_over_rho_coeff); // actio = - reactio (but twice! nothing changes)

        // compute and add the correct densityDot_j
        LINALG::Assemble(*densityDotn, densityDotn_j, particle_j.gid, particle_j.owner);
        // compute and add the correct specEnthalpyDot_j
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_j, particle_j.gid, particle_j.owner);
        // compute and add the correct accn_j. action-reaction, there is += instead of -
        LINALG::Assemble(*accn, accn_j, particle_j.lm, particle_j.owner);
        // compute and add the correct gradT_over_rho_i
        //LINALG::Assemble(*gradT_over_rho, gradT_over_rho_j, particle_j.lm, particle_j.owner);

      }
    }
  }
}

/*----------------------------------------------------------------------*
 | compute interaction with neighboring walls              katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringWallMeshFreeInteraction(
    const int lidNodeCol_i,
    const std::set<DRT::Element*>& neighboring_walls,
    const Teuchos::RCP<DRT::Discretization>& walldiscret,
    const Teuchos::RCP<const Epetra_Vector>& walldisn,
    const Teuchos::RCP<const Epetra_Vector>& wallveln,
    const Teuchos::RCP<Epetra_Vector>& accn)
{

  // ids i
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];

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
  std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(discret_->lColNode(lidNodeCol_i))->Get_history_wall();
  if(history_wall.size() > 3)
    dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

  const double p_Rho2_i = particle_i.pressure/std::pow(particle_i.density,2);

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    WallInteractionPoint wallcontact = surfaces[s];

    // compute distance and relative velocities
    LINALG::Matrix<3,1> rRel;
    double WFDerivative = 0;
    rRel.Update(1.0, particle_i.dis, -1.0, wallcontact.point);
    const double rRelNorm2 = rRel.Norm2();
    LINALG::Matrix<3,1> vRel(particle_i.vel);
    const double rRelDotVrel = rRel.Dot(vRel);

    if(rRelNorm2 == 0.0)
      dserror("particle center and wall are lying in the same place -> bad initialization?");

    // compute the proper weight function derivative
    switch (weightFunctionType_)
    {
    case INPAR::PARTICLE::CubicBspline :
    {
      WFDerivative = PARTICLE::WeightFunction_CubicBspline::DerivativeWeight(rRelNorm2, particle_i.radius);
      break;
    }
    }

    if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
    {
      wallMeshFreeData_.pressure = particle_i.pressure;
      wallMeshFreeData_.density = particle_i.density;
      wallMeshFreeData_.mass = particle_i.mass;
    }

    // p_i/\rho_i^2 + p_j/\rho_j^2
    const double gradP_Rho2 = (p_Rho2_i + wallMeshFreeData_.pressure/(wallMeshFreeData_.density * wallMeshFreeData_.density));
    // WFDerivative/(\rho_i \rho_j)
    const double WFDerivative_Rho2 = WFDerivative / (particle_i.density * wallMeshFreeData_.density);

    LINALG::Matrix<3,1> momentum_i;
    // compute the pressure term
    momentum_i.Update(- gradP_Rho2 * WFDerivative,rRel);
    // compute the diffusion term
    momentum_i.Update(diffusionCoeff_ * WFDerivative_Rho2,vRel,1.0);
    // compute the convection term
    momentum_i.Update(convectionCoeff_ * WFDerivative_Rho2 * rRelDotVrel / (rRelNorm2 * rRelNorm2),rRel,1.0);

    // mass scalings
    LINALG::Matrix<3,1> accn_i;
    accn_i.Update(wallMeshFreeData_.mass,momentum_i);

    LINALG::Assemble(*accn, accn_i, particle_i.lm, particle_i.owner);
  }
}


/*--------------------------------------------------------------------------*
 | calculate interaction with neighboring heat sources         katta 10/16  |
 *--------------------------------------------------------------------------*/
void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringHeatSourcesContact(
  const int lidNodeCol_i,
  const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources,
  const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn)
{
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];

  double specEnthalpyDot_i = 0.0;
  std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less>::const_iterator hs;
  for(hs = neighboring_heatSources->begin(); hs != neighboring_heatSources->end();  ++hs)
  {
    if ((*hs)->minVerZone_[0]<=particle_i.dis(0) &&
        (*hs)->minVerZone_[1]<=particle_i.dis(1) &&
        (*hs)->minVerZone_[2]<=particle_i.dis(2) &&
        (*hs)->maxVerZone_[0]>=particle_i.dis(0) &&
        (*hs)->maxVerZone_[1]>=particle_i.dis(1) &&
        (*hs)->maxVerZone_[2]>=particle_i.dis(2))
    {
      specEnthalpyDot_i += ((*hs)->QDot_)/particle_i.density;
    }
  }

  LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDot_i, particle_i.gid, particle_i.owner);
}

/*----------------------------------------------------------------------*
 | calc neighbouring particleMeshFree interactions         katta 10/16  |
 *----------------------------------------------------------------------*/
/// --- second round -> unactivated ----------------------------------- ///
/*
void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringParticleMeshFreeInteraction_SecondRound(
  const int lidNodeCol_i,
  const std::list<DRT::Node*> neighboring_particles,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn)
{
  // ids i
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];

  // extract the interesting material properties
  const double thermalConductivity = particle_algorithm_->ExtParticleMat()->thermalConductivity_;

  // loop over the neighbouring particles
  // avoid to recompute forces (self-interaction are not allowed in this part of the code)
  for(std::list<DRT::Node*>::const_iterator jj=neighboring_particles.begin(); jj!=neighboring_particles.end(); ++jj)
  {
    const ParticleMeshFreeData &particle_j = particleMeshFreeData_[(*jj)->LID()];

    // evaluate contact only once in case we own particle j and the radii match. Otherwise compute everything,
    // the assemble method does not write in case the particle is a ghost.
    // another check on the radii is performed (later in the code) to handle the case where we want to write, we can write but accelerations are not actio<->reactio
    if(particle_i.gid < particle_j.gid || particle_i.radius != particle_j.radius || particle_j.owner != myrank_)
    {

      // compute distance and relative velocities
      LINALG::Matrix<3,1> rRel, gradT_over_rhoRel, WFGrad;
      rRel.Update(1.0, particle_i.dis, -1.0, particle_j.dis);
      gradT_over_rhoRel.Update(1.0, particle_i.gradT_over_rho, 1.0, particle_j.gradT_over_rho); // divergence = sum

      // compute the proper weight function gradient
      switch (weightFunctionType_)
      {
      case INPAR::PARTICLE::CubicBspline :
      {
        WFGrad = PARTICLE::WeightFunction_CubicBspline::GradientWeight(rRel, particle_i.radius);
        break;
      }
      }
      // compute WFGradDotGradT_over_rhoRel
      const double WFGradDotGradT_over_rhoRel = WFGrad.Dot(gradT_over_rhoRel);

      double specEnthalpyDotn_i = thermalConductivity * particle_j.mass * WFGradDotGradT_over_rhoRel / particle_j.density;

      // compute and add the correct densityDot_i
      LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_i, particle_i.gid, particle_i.owner);

      // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
      // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
      if(particle_i.radius == particle_j.radius)
      {
        double specEnthalpyDotn_j = - thermalConductivity * particle_i.mass * WFGradDotGradT_over_rhoRel / particle_i.density; // actio = - reactio

        // compute and add the correct densityDot_j
        // beware! we assumed that WFGrad_ij = - WFGrad_ji and vRel_ij = - vRel_ji
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_j, particle_j.gid, particle_j.owner);
      }
    }
  }
}
*/
///--------------------------------------------------------------------///
