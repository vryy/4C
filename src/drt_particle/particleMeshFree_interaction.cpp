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

  // vectors to collect the inputs of the second round

  Teuchos::RCP<Epetra_Vector> colorFieldGradientn = Teuchos::null;
  // list of the list of all particles in the neighborhood of currparticle and bins for multiple iterations
  std::list<DRT::Element*> binList;
  std::list<std::list<DRT::Node*> > neighboring_particlesList;
  // create the vector only in case it is needed
  if (trg_secondRound_)
  {
    colorFieldGradientn = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }

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

    if (trg_secondRound_)
    {
      // store for future use
      binList.push_back(currentBin);
      neighboring_particlesList.push_back(neighboring_particles);
    }

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing
      DRT::Node* particle_i = currentBinParticles[i];
      const int lidNodeCol_i = particle_i->LID();

      // compute neighbouring heat sources
      CalcNeighboringHeatSourcesContact(lidNodeCol_i, neighboring_heatSources, specEnthalpyDotn);

      // compute contact with neighboring walls
      CalcNeighboringWallMeshFreeInteraction(lidNodeCol_i, neighboring_walls, walldiscret, walldisn, wallveln, accn, densityDotn);

      // compute interactions with neighboring particles
      CalcNeighboringParticleMeshFreeInteraction(lidNodeCol_i, neighboring_particles, accn, densityDotn, specEnthalpyDotn, colorFieldGradientn);
    }
  }

  // second round
  if (trg_secondRound_)
  {
    // reorder the interesting vectors
    SetStateVectors_SecondRound(colorFieldGradientn);

    std::list<std::list<DRT::Node*> >::const_iterator currNeighboring_particles = neighboring_particlesList.begin();
    for (std::list<DRT::Element*>::const_iterator currBin = binList.begin(); currBin != binList.end(); ++currBin)
    {
      // extract the pointer to the particles
      DRT::Node** currentBinParticles = (*currBin)->Nodes();

      // loop over all particles in CurrentBin
      for(int i=0; i<(*currBin)->NumNode(); ++i)
      {
        // determine the particle we are analizing
        DRT::Node* particle_i = currentBinParticles[i];
        const int lidNodeCol_i = particle_i->LID();

        // compute interactions with neighboring particles - second round
        CalcNeighboringParticleMeshFreeInteraction_SecondRound(lidNodeCol_i, (*currNeighboring_particles), accn);
      }
      ++currNeighboring_particles;
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
  Teuchos::RCP<Epetra_Vector> densityDotn,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
  Teuchos::RCP<Epetra_Vector> colorFieldGradientn)
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
      LINALG::Matrix<3,1> rRel, rRelVersor;
      rRel.Update(1.0, particle_i.dis, -1.0, particle_j.dis); // inward vector
      const double rRelNorm2 = rRel.Norm2();

      // skip in case particles are too apart
      if (rRelNorm2>=particle_i.radius)
      {
        continue;
      }

      // --- general quantities --- //

      rRelVersor.Update(1.0/rRelNorm2,rRel);

      LINALG::Matrix<3,1> vRel;
      vRel.Update(1.0, particle_i.vel, -1.0, particle_j.vel);

      const double rRelVersorDotVrel = rRelVersor.Dot(vRel);

      // compute the proper weight function gradient
      const double weightDerivative = weightDerivative_(rRelNorm2, particle_i.radius);

      // correction parameter for adhesion-cohesion mechanics
      const double densityCorrectiveTerm = 2 * restDensity_ / (particle_i.density + particle_j.density);

      // --- density --- //

      const double WFGradDotVrel = weightDerivative * rRelVersorDotVrel;

      // mass scaling and assembling
      double densityDotn_i = particle_j.mass * WFGradDotVrel;
      LINALG::Assemble(*densityDotn, densityDotn_i, particle_i.gid, particle_i.owner);

      // --- acceleration --- //

      LINALG::Matrix<3,1> momentum_i;

      // compute the pressure term
      const double rhoSquare_j = std::pow(particle_j.density,2);
      const double gradP_Rho2 = (p_Rho2_i + particle_j.pressure/rhoSquare_j);  // p_i/\rho_i^2 + p_j/\rho_j^2
      momentum_i.Update(- gradP_Rho2 * weightDerivative,rRelVersor);

      // compute the diffusion term
      const double weightDerivative_Rho2 = weightDerivative / (particle_i.density * particle_j.density);
      momentum_i.Update(diffusionCoeff_ * weightDerivative_Rho2/rRelNorm2,vRel,1.0);

      // compute the convection term
      momentum_i.Update(convectionCoeff_ * weightDerivative_Rho2 * rRelVersorDotVrel / rRelNorm2,rRelVersor,1.0);

      // compute the cohesion term
      const double cohesionWeight = PARTICLE::SurfaceTensionInteractions::Cohesion(rRelNorm2, particle_i.radius);
      momentum_i.Update( - densityCorrectiveTerm * surfaceVoidTension_ * cohesionWeight,rRelVersor,1.0);

      // mass scaling and assembling
      LINALG::Matrix<3,1> accn_i;
      accn_i.Update(particle_j.mass, momentum_i);
      LINALG::Assemble(*accn, accn_i, particle_i.lm, particle_i.owner);

      // --- specific enthalpy --- //

      const double divT_rho_i = 2 * particle_algorithm_->ExtParticleMat()->thermalConductivity_ * (weightDerivative_Rho2/rRelNorm2) * (particle_i.temperature - particle_j.temperature);

      // mass scaling and assembling
      double specEnthalpyDotn_i = particle_j.mass * divT_rho_i;
      LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_i, particle_i.gid, particle_i.owner);

      // --- color field --- //

      LINALG::Matrix<3,1> colorFieldGradient_mass_i;
      if (trg_secondRound_)
      {
        colorFieldGradient_mass_i.Update(particle_i.radius * weightDerivative / particle_j.density,rRelVersor);

        // mass scaling and assembling
        LINALG::Matrix<3,1> colorFieldGradientn_i;
        colorFieldGradientn_i.Update(particle_j.mass, colorFieldGradient_mass_i);
        LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_i, particle_i.lm, particle_i.owner);
      }

      // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
      // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
      if(particle_i.radius == particle_j.radius)
      {
        // --- density --- //

        // mass scaling and assembling
        double densityDotn_j = particle_i.mass * WFGradDotVrel; // actio = - reactio (but twice! nothing changes)
        LINALG::Assemble(*densityDotn, densityDotn_j, particle_j.gid, particle_j.owner);

        // --- acceleration --- //

        // mass scaling and assembling
        LINALG::Matrix<3,1> accn_j;
        accn_j.Update(- particle_i.mass,momentum_i); // actio = - reactio -> momentum_j = - momentum_i
        LINALG::Assemble(*accn, accn_j, particle_j.lm, particle_j.owner); // action-reaction, there is += instead of -

        // --- specific enthalpy --- //

        // mass scaling and assembling
        double specEnthalpyDotn_j = - particle_i.mass * divT_rho_i; // actio = - reactio -> divT_rho_j = - divT_rho_i
        LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn_j, particle_j.gid, particle_j.owner);

        // --- color field --- //

        if (trg_secondRound_)
        {
          // mass scaling and assembling
          LINALG::Matrix<3,1> colorFieldGradientn_j;
          colorFieldGradientn_j.Update(- particle_i.mass, colorFieldGradient_mass_i); // actio = - reactio -> colorFieldGradientn_j = colorFieldGradientn_i
          LINALG::Assemble(*colorFieldGradientn, colorFieldGradientn_j, particle_j.lm, particle_j.owner);
        }
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
    const Teuchos::RCP<Epetra_Vector>& accn,
    Teuchos::RCP<Epetra_Vector> densityDotn)
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

  const double rhoSquare_i = std::pow(particle_i.density,2);
  const double p_Rho2_i = particle_i.pressure/rhoSquare_i;

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    WallInteractionPoint wallcontact = surfaces[s];

    // compute distance and relative velocities
    LINALG::Matrix<3,1> rRel, rRelVersor;
    rRel.Update(1.0, particle_i.dis, -1.0, wallcontact.point);
    const double rRelNorm2 = rRel.Norm2();

    // skip in case particles are too apart
    if (rRelNorm2>=particle_i.radius)
    {
      continue;
    }

    rRelVersor.Update(1/rRelNorm2,rRel);

    LINALG::Matrix<3,1> vRel(particle_i.vel);

    if (wallInteractionType_ == INPAR::PARTICLE::Mirror)
    {
      wallMeshFreeData_.pressure = particle_i.pressure;
      wallMeshFreeData_.density = particle_i.density;
      wallMeshFreeData_.mass = particle_i.mass;
      // mirrored velocity
      vRel.Scale(2);
    }

    const double rRelVersorDotVrel = rRelVersor.Dot(vRel);

    if(rRelNorm2 == 0.0)
    {
      dserror("particle center and wall are lying in the same place -> bad initialization?");
    }

    // --- general quantities --- //

    // compute the proper weight function derivative
    const double weightDerivative = weightDerivative_(rRelNorm2, particle_i.radius);

    // --- density --- //

    const double WFGradDotVrel = weightDerivative * rRelVersorDotVrel; // compute WFGradDotVrel

    // mass scaling and assembling
    double densityDotn_i = wallMeshFreeData_.mass * WFGradDotVrel;
    LINALG::Assemble(*densityDotn, densityDotn_i, particle_i.gid, particle_i.owner);

    // --- acceleration --- //

    LINALG::Matrix<3,1> momentum_i;

    // compute the pressure term
    const double gradP_Rho2 = (p_Rho2_i + wallMeshFreeData_.pressure/(wallMeshFreeData_.density * wallMeshFreeData_.density));  // p_i/\rho_i^2 + p_j/\rho_j^2
    momentum_i.Update(- gradP_Rho2 * weightDerivative,rRelVersor);

    // compute the diffusion term
    const double weightDerivative_Rho2 = weightDerivative / (particle_i.density * wallMeshFreeData_.density); // WFDerivative/(\rho_i \rho_j)
    momentum_i.Update(diffusionCoeff_ * weightDerivative_Rho2/rRelNorm2,vRel,1.0);

    // compute the convection term
    momentum_i.Update(convectionCoeff_ * weightDerivative_Rho2 * rRelVersorDotVrel / rRelNorm2,rRelVersor,1.0);

    // compute the adhesion term
    const double adhesionWeight = PARTICLE::SurfaceTensionInteractions::Adhesion(rRelNorm2, particle_i.radius);
    momentum_i.Update(- surfaceWallTension_ * adhesionWeight,rRelVersor,1.0);

    // mass scaling and assembling
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

void PARTICLE::ParticleMeshFreeInteractionHandler::CalcNeighboringParticleMeshFreeInteraction_SecondRound(
  const int lidNodeCol_i,
  const std::list<DRT::Node*> neighboring_particles,
  Teuchos::RCP<Epetra_Vector> accn)
{
  // ids i
  const ParticleMeshFreeData &particle_i = particleMeshFreeData_[lidNodeCol_i];


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
      // --- general quantities --- //

      // correction parameter for adhesion-cohesion mechanics
      const double densityCorrectiveTerm = 2 * restDensity_ / (particle_i.density + particle_j.density);

      // --- surface tension - colorFieldGradient --- //

      // rescaling and assembpling
      LINALG::Matrix<3,1> accn_i(particle_j.colorFieldGradient);
      accn_i.Update(1.0, particle_i.colorFieldGradient, -1.0);
      accn_i.Scale(- densityCorrectiveTerm * surfaceVoidTension_);
      LINALG::Assemble(*accn, accn_i, particle_i.lm, particle_i.owner);

      // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
      // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
      if(particle_i.radius == particle_j.radius)
      {
        // --- surface tension - colorFieldGradient --- //

        LINALG::Matrix<3,1> accn_j;
        accn_j.Update(-1.0, accn_i); // actio = - reactio
        LINALG::Assemble(*accn, accn_j, particle_j.lm, particle_j.owner);
      }
    }
  }
}

