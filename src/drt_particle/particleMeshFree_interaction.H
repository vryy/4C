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

    //! set state vectors (col type)
    virtual void SetStateVectors();

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual void EvaluateParticleMeshFreeInteractions(
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<Epetra_Vector> densityDotn);

    virtual void CalcNeighboringParticleMeshFreeInteraction(
      const int lidNodeCol_i,
      const std::list<DRT::Node*> neighboring_particles,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<Epetra_Vector> densityDotn);

    /// calculate interaction with neighboring walls
    virtual void CalcNeighboringWallMeshFreeInteraction(
      const int lidNodeCol_i,
      const std::set<DRT::Element*>& neighboring_walls,
      const Teuchos::RCP<DRT::Discretization>& walldiscret,
      const Teuchos::RCP<const Epetra_Vector>& walldisn,
      const Teuchos::RCP<const Epetra_Vector>& wallveln,
      const Teuchos::RCP<Epetra_Vector>& accn);

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

    /// temporary storage of particleMeshFree data
    std::vector<ParticleMeshFreeData> particleMeshFreeData_;

    /// weight function type
    const INPAR::PARTICLE::WeightFunction weightFunctionType_;

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    //@}

    //! @name wall interaction parameters
    //@{


    /// pressure divergence of the wall
    double wallInteractionPressureDivergence_;

    ///fake mass of the virtual particles on the wall
    double wallInteractionFakeMass_;

    //@}

  };  // class ParticleCollisionHandler

  //! @name Container structs for faster access
  //@{

  /// struct to store particleMeshFree data for faster access
  struct ParticleMeshFreeData {
    LINALG::Matrix<3,1> dis;
    LINALG::Matrix<3,1> vel;
    double radius;
    double density;
    double mass;
    double pressure;
    int gid;
    std::vector<int> lm;
    int owner;
  };

  /// struct to store wall contact point and additional information about corresponding wall ele
  struct WallInteractionPoint {
   int elemId;
   LINALG::Matrix<3,1> point;
   double penetration;
   std::map<int,LINALG::Matrix<3,1> > nodalCoordinates;
   std::vector<int> lm;
   std::vector<int> lmowner;
  };

  //@}

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
