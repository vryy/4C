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

    //@}


  protected:

    //! @name weight functions
    //@{

    /// compute the weight function gradient (cubicBspline type)
    virtual void GradientWeightFunction_CubicBSpline(LINALG::Matrix<3,1> &WFGrad, const double &radius);

    //@}

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

  //@}

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
