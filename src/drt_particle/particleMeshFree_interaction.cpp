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
#include "particleMeshFree_weightFunction.H"
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
  myrank_(discret->Comm().MyPID())
{
// checks
if (particle_algorithm_->ExtParticleMat() == NULL)
  dserror("extParticleMat_ is empty");

// extract the timint object
const Teuchos::RCP<ADAPTER::Particle> adapterParticle = particle_algorithm_->AdapterParticle();
// fill the col vectors with the proper particle data
SetStateVectors(
  adapterParticle->Dispnp(),
  adapterParticle->Velnp(),
  adapterParticle->Radiusnp(),
  adapterParticle->Densitynp(),
  adapterParticle->Pressure());

// set the proper weightFunction object
const INPAR::PARTICLE::WeightFunction weightFunctionType = DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION");
switch (weightFunctionType)
{
  case INPAR::PARTICLE::CubicBspline :
  {
    weightFunction_ = Teuchos::rcp(new PARTICLE::WeightFunction_CubicBspline());
    break;
  }
}

}

/*----------------------------------------------------------------------*
 | set colVectors                                          katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::SetStateVectors(
  Teuchos::RCP<const Epetra_Vector> disn,
  Teuchos::RCP<const Epetra_Vector> veln,
  Teuchos::RCP<const Epetra_Vector> radiusn,
  Teuchos::RCP<const Epetra_Vector> densityn,
  Teuchos::RCP<const Epetra_Vector> pressure)
{
// checks
if (disn == Teuchos::null     ||
    veln == Teuchos::null     ||
    densityn == Teuchos::null ||
    pressure == Teuchos::null)
  dserror("one or more state vectors are empty");

/// miraculous transformation into column vectors... ///

// dof based vectors
Teuchos::RCP<Epetra_Vector> disnCol_ = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*disn,*disnCol_);
Teuchos::RCP<Epetra_Vector> velnCol_ = LINALG::CreateVector(*discret_->DofColMap(),false);
LINALG::Export(*veln,*velnCol_);
// node based vectors
Teuchos::RCP<Epetra_Vector> radiusnCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*radiusn,*radiusnCol);
Teuchos::RCP<Epetra_Vector> densitynCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*densityn,*densitynCol);
Teuchos::RCP<Epetra_Vector> pressureCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
LINALG::Export(*pressure,*pressureCol);
}

/*----------------------------------------------------------------------*
 | evaluate interactions                                   katta 10/16  |
 *----------------------------------------------------------------------*/

void PARTICLE::ParticleMeshFreeInteractionHandler::EvaluateParticleMeshFreeInteractions(
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<Epetra_Vector> densityDotn)
{
  std::cout << "Work in progress!!!\n";
}



