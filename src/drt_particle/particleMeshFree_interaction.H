/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.H

\brief Particle-MeshFree interaction handling

\level 3

\maintainer Alessandro Cattabiani
*/

/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_
#define SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_
/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle.H"
#include "particle_node.H"
#include "../drt_lib/drt_element.H"
#include <Epetra_FEVector.h>


// forward declarations
namespace DRT
{
  class Element;
  class Node;
}
namespace PARTICLE
{
  class Algorithm;
  class WeightFunctionBase;
  struct WallContactPoint;
  struct ParticleMeshFreeData;
}

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{
  /*====================================================================*/
  /*!
   * \brief Interaction handler for particle-meshFree interactions
   *
   * \author cattabiani
   * \date 10/16
   */
  class ParticleMeshFreeInteractionHandler
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    ParticleMeshFreeInteractionHandler(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
      const Teuchos::ParameterList& particledynparams);

    //! Destructor
    virtual ~ParticleMeshFreeInteractionHandler(){ ; }

    //@}

    //! @name Actions
    //@{

    //! set col vectors
    virtual void SetStateVectors(
      Teuchos::RCP<const Epetra_Vector> disn,
      Teuchos::RCP<const Epetra_Vector> veln,
      Teuchos::RCP<const Epetra_Vector> radiusn,
      Teuchos::RCP<const Epetra_Vector> densityn,
      Teuchos::RCP<const Epetra_Vector> pressure);

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual void EvaluateParticleMeshFreeInteractions(
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<Epetra_Vector> densityDotn);

    //@}


  protected:

    //! @name helper classes
    //@{

    ///< particle discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    ///< particle algorithm that is connected to the time integration
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //@}

    //! @name components
    //@{

    /// ID of actual processor in parallel
    int myrank_;

    /// Weight function
    Teuchos::RCP<PARTICLE::WeightFunctionBase> weightFunction_;

    //@}

    //! @name col vectors
    //@{

    Teuchos::RCP<Epetra_Vector> disnCol_;
    Teuchos::RCP<Epetra_Vector> velnCol_;
    Teuchos::RCP<Epetra_Vector> radiusnCol_;
    Teuchos::RCP<Epetra_Vector> densitynCol_;
    Teuchos::RCP<Epetra_Vector> pressureCol_;

    //@}

  };  // class ParticleCollisionHandler

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
