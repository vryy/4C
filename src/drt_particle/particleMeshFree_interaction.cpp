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
  weightFunctionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION"))
{
// checks
if (particle_algorithm_->ExtParticleMat() == NULL)
  dserror("extParticleMat_ is empty");


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

  std::cout << "wall contact is not implemented!\n";
  //const bool havepbc = particle_algorithm_->HavePBCs();

  // store bins, which have already been examined
  std::set<int> examinedbins;

  // list of all particles in the neighborhood of currparticle
  std::list<DRT::Node*> neighboring_particles;

  // loop over the particles (no superpositions)
  const int numrowparticles = discret_->NodeRowMap()->NumMyElements();

  for(int i_skippingLoop=0; i_skippingLoop<numrowparticles; ++i_skippingLoop)
  {
    // extract the particle
    DRT::Node *currparticle = discret_->lRowNode(i_skippingLoop);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    // check I own the bin
    assert(currentBin->Owner() == myrank_);

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
      //CalcNeighboringWallsContact(particle_i, data_i, neighboring_walls, dt,
      //    walldiscret, walldisn, wallveln, f_contact, m_contact, f_structure);

      // compute interactions with neighboring particles
      CalcNeighboringParticleMeshFreeInteraction(lidNodeCol_i,neighboring_particles,accn,densityDotn);

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
  const int lidNodeRow_i = densityDotn->Map().LID(particle_i.gid);
  const double p_over_rhoSquare_i = particle_i.pressure/std::pow(particle_i.density,2);

  // self-interaction
  // densityDot -> the weightFunction gradient is null (or ill-posed in case of a strange weightFunction)
  // acc -> the weightFunction gradient is null (or ill-posed in case of a strange weightFunction)

  // loop over the neighbouring particles
  // avoid to recompute forces (self-interaction are not allowed in this part of the code)
  for(std::list<DRT::Node*>::const_iterator jj=neighboring_particles.begin(); jj!=neighboring_particles.end(); ++jj)
  {
    const ParticleMeshFreeData &particle_j = particleMeshFreeData_[(*jj)->LID()];

    // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
    // if you pass this checkpoint you can compute the effect of particle j on particle i
    if(particle_i.gid >= particle_j.gid && particle_i.radius == particle_j.radius && particle_j.owner == myrank_)
      continue;

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

    // compute and add the correct densityDot_i
    (*densityDotn)[lidNodeRow_i] += particle_j.mass * WFGradDotVrel;

    // p_i/\rho_i^2 + p_j/\rho_j^2
    const double divergenceCoeff = (p_over_rhoSquare_i + particle_j.pressure/std::pow(particle_j.density,2));
    // compute and add the correct accn_i
    for (int dim = 0; dim<3; ++dim)
      (*accn)[3*lidNodeRow_i + dim] -= particle_j.mass * divergenceCoeff * WFGrad(dim);

    // evaluate contact only once if possible! (the logic table has been double-checked, trust me)
    // if you are here and you pass this checkpoint you can compute the effect of particle i on particle j
    if(!(particle_i.gid < particle_j.gid && particle_i.radius == particle_j.radius && particle_j.owner == myrank_))
      continue;

    const int lidNodeRow_j = densityDotn->Map().LID(particle_j.gid);

    // compute and add the correct densityDot_j
    // beware! we assumed that WFGrad_ij = - WFGrad_ji and vRel_ij = - vRel_ji
    (*densityDotn)[lidNodeRow_j] += particle_i.mass * WFGradDotVrel;
    // compute and add the correct accn_j. action-reaction, there is += instead of -
    for (int dim = 0; dim<3; ++dim)
      (*accn)[3*lidNodeRow_j + dim] += particle_i.mass * divergenceCoeff * WFGrad(dim);
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


