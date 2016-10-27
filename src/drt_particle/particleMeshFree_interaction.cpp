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
#include "particle_algorithm.H"
#include "particle_timint_centrdiff.H"
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
switch(wallInteractionType_)
{
case INPAR::PARTICLE::InitParticle :
  {
    const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
    wallInteractionPressureDivergence_ = std::pow(extParticleMat->SpeedOfSoundS(),2) / extParticleMat->initDensity_;
    wallInteractionFakeMass_ = extParticleMat->initDensity_ * (4.0 / 3.0) * M_PI * std::pow(extParticleMat->initRadius_,3);
    break;
  }
case INPAR::PARTICLE::Mirror :
  {
    wallInteractionPressureDivergence_ = -1;
    wallInteractionFakeMass_ = -1;
    break;
  }
case INPAR::PARTICLE::Custom :
  {
    const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
    wallInteractionPressureDivergence_ = particleparams.get<double>("WALL_INTERACTION_PRESSDIV");
    wallInteractionFakeMass_ = particleparams.get<double>("WALL_INTERACTION_FAKEMASS");
  }
}

// other checks
if (wallInteractionType_ != INPAR::PARTICLE::Mirror)
{
  if (wallInteractionPressureDivergence_ < 0)
    dserror("the value of wallInteractionPressureDivergence_ is unacceptable");
  if (wallInteractionFakeMass_ < 0)
    dserror("the value of wallInteractionFakeMass_ is unacceptable");
}
// fill the col vectors with the proper particle data
SetStateVectors();

}

/*----------------------------------------------------------------------*
 | set colVectors                                          katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVectors()
{
  // extract the timint object
  const Teuchos::RCP<ADAPTER::Particle> adapterParticle = particle_algorithm_->AdapterParticle();

  Teuchos::RCP<const Epetra_Vector> disn = adapterParticle->Dispnp();
  Teuchos::RCP<const Epetra_Vector> veln = adapterParticle->Velnp();
  Teuchos::RCP<const Epetra_Vector> radiusn = adapterParticle->Radiusnp();
  Teuchos::RCP<const Epetra_Vector> densityn = adapterParticle->Densitynp();
  Teuchos::RCP<const Epetra_Vector> mass = adapterParticle->Mass();
  Teuchos::RCP<const Epetra_Vector> pressure = adapterParticle->Pressure();

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
Teuchos::RCP<Epetra_Vector> massCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*mass,*massCol);
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
  data.mass = (*massCol)[lid];
  data.pressure = (*pressureCol)[lid];

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
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::EvaluateParticleMeshFreeInteractions(
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<Epetra_Vector> densityDotn)
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

  // list of all particles in the neighborhood of currparticle
  std::list<DRT::Node*> neighboring_particles;

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

    // remove current content but keep memory
    neighboring_particles.clear();

    // list of walls that border on the CurrentBin
    std::set<DRT::Element*> neighboring_walls;

    particle_algorithm_->GetNeighbouringParticlesAndWalls(binId, neighboring_particles, neighboring_walls);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing
      DRT::Node* particle_i = currentBinParticles[i];
      const int lidNodeCol_i = particle_i->LID();

      // compute contact with neighboring walls
      CalcNeighboringWallMeshFreeInteraction(lidNodeCol_i, neighboring_walls, walldiscret, walldisn, wallveln, accn);

      // compute interactions with neighboring particles
      CalcNeighboringParticleMeshFreeInteraction(lidNodeCol_i, neighboring_particles, accn, densityDotn);

    }
  }
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
  Teuchos::RCP<Epetra_Vector> densityDotn)
{
  // ids i
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];
  const double p_over_rhoSquare_i = particle_i.pressure/std::pow(particle_i.density,2);

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
      LINALG::Matrix<3,1> WFGrad, vRel;
      WFGrad.Update(1.0, particle_i.dis, -1.0, particle_j.dis);
      vRel.Update(1.0, particle_i.vel, -1.0, particle_j.vel);

      // compute the proper weight function gradient
      switch (weightFunctionType_)
      {
      case INPAR::PARTICLE::CubicBspline :
      {
        GradientWeightFunction_CubicBSpline(WFGrad, particle_i.radius);
        break;
      }
      }

      // compute WFGradVrel
      const double WFGradDotVrel = WFGrad.Dot(vRel);
      // p_i/\rho_i^2 + p_j/\rho_j^2
      const double divergenceCoeff = (p_over_rhoSquare_i + particle_j.pressure/std::pow(particle_j.density,2));

      static std::vector<int> nodeOwner_i(1);
      nodeOwner_i[0] = particle_i.owner;
      static std::vector<int> nodeGid_i(1);
      nodeGid_i[0] = particle_i.gid;
      static Epetra_SerialDenseVector densityDotn_i(1);
      densityDotn_i[0] = particle_j.mass * WFGradDotVrel;

      static Epetra_SerialDenseVector accn_i(3);
      static std::vector<int> lmowner_i(3);
      for(unsigned dim=0; dim<3; ++dim)
      {
        accn_i[dim] = - particle_j.mass * divergenceCoeff * WFGrad(dim);
        lmowner_i[dim] = particle_i.owner;
      }

      // compute and add the correct densityDot_i
      LINALG::Assemble(*densityDotn, densityDotn_i, nodeGid_i, nodeOwner_i);
      // compute and add the correct accn_i
      LINALG::Assemble(*accn, accn_i, particle_i.lm, lmowner_i);

      // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
      // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
      if(particle_i.radius == particle_j.radius)
      {
        static std::vector<int> nodeOwner_j(1);
        nodeOwner_j[0] = particle_j.owner;
        static std::vector<int> nodeGid_j(1);
        nodeGid_j[0] = particle_j.gid;
        static Epetra_SerialDenseVector densityDotn_j(1);
        densityDotn_j[0] = particle_i.mass * WFGradDotVrel;

        static Epetra_SerialDenseVector accn_j(3);
        static std::vector<int> lmowner_j(3);
        for(unsigned dim=0; dim<3; ++dim)
        {
          accn_j[dim] = particle_i.mass * divergenceCoeff * WFGrad(dim); //there is no - because: actio = - reactio
          lmowner_j[dim] = particle_j.owner;
        }

        // compute and add the correct densityDot_j
        // beware! we assumed that WFGrad_ij = - WFGrad_ji and vRel_ij = - vRel_ji
        LINALG::Assemble(*densityDotn, densityDotn_j, nodeGid_j, nodeOwner_j);
        // compute and add the correct accn_j. action-reaction, there is += instead of -
        LINALG::Assemble(*accn, accn_j, particle_j.lm, lmowner_j);
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

  const double p_over_rhoSquare_i = particle_i.pressure/std::pow(particle_i.density,2);

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    WallInteractionPoint wallcontact = surfaces[s];

    // distance-vector
    static LINALG::Matrix<3,1> WFGrad;
    WFGrad.Update(1.0, particle_i.dis, -1.0, wallcontact.point);

    if(WFGrad.Norm2() == 0.0)
      dserror("particle center and wall are lying in the same place -> bad initialization?");

    // compute the proper weight function gradient
    switch (weightFunctionType_)
    {
    case INPAR::PARTICLE::CubicBspline :
    {
      GradientWeightFunction_CubicBSpline(WFGrad, particle_i.radius);
      break;
    }
    }

    if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
    {
      wallInteractionPressureDivergence_ = p_over_rhoSquare_i;
      wallInteractionFakeMass_ = particle_i.mass;
    }

    // p_i/\rho_i^2 + p_j/\rho_j^2
    const double divergenceCoeff = (p_over_rhoSquare_i + wallInteractionPressureDivergence_);

    static Epetra_SerialDenseVector accn_i(3);
    static std::vector<int> lmowner_i(3);
    for(unsigned dim=0; dim<3; ++dim)
    {
      accn_i[dim] = - wallInteractionFakeMass_ * divergenceCoeff * WFGrad(dim);
      lmowner_i[dim] = particle_i.owner;
    }

    // compute and add the correct accn_i
    LINALG::Assemble(*accn, accn_i, particle_i.lm, lmowner_i);
  }
}

/*-------------------------------------------------------------------------------*
 | compute the weight function gradient (cubicBspline type)         katta 10/16  |
 *-------------------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::GradientWeightFunction_CubicBSpline(LINALG::Matrix<3,1> &WFGrad, const double &radius)
{
  // safety checks
  assert(radius > 0);

  const double norm_const = 1.5 * M_1_PI;
  const double WFGradNorm = WFGrad.Norm2();

  // solving the particular case in which two particles perfectly overlap
  if (WFGradNorm <= 1e-16)
  {
    dserror("Warning! particles are overlapping! Right now, it is not allowed");
    std::cout << "Warning! particles are overlapping!\n";
    WFGrad.PutScalar(0);
    return;
  }

  const double resizer_temp = 2 / radius;
  const double norm_dist_rel = resizer_temp * WFGradNorm;
  const double resizer = resizer_temp / WFGradNorm;

  if (norm_dist_rel< 1)
    WFGrad.Scale(norm_const * resizer * (1.5 * std::pow(norm_dist_rel,2) - 2 * norm_dist_rel));
  else if (norm_dist_rel< 2)
    WFGrad.Scale(- 0.5 * norm_const * resizer * std::pow(norm_dist_rel - 2,2));
  else
    WFGrad.PutScalar(0);

}
