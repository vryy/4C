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

  /// struct to store wallMeshFree data for faster access
  struct WallMeshFreeData {
    double density;
    double mass;
    double pressure;
  };

  class HeatSource;
}

namespace BINSTRATEGY
{
  class Less;
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
    virtual void SetStateVectors(
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> pressure);
    /// --- second round -> unactivated ----------------------------------- ///
    /*
    //! set state vectors (col type) - second round
    virtual void SetStateVectors_SecondRound(
        Teuchos::RCP<const Epetra_Vector> gradT_over_rho);

    /// calculate interaction with neighboring particles - second round
    virtual void CalcNeighboringParticleMeshFreeInteraction_SecondRound(
      const int lidNodeCol_i,
      const std::list<DRT::Node*> neighboring_particles,
      Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);
    */
    ///--------------------------------------------------------------------///

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual void EvaluateParticleMeshFreeInteractions(
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<Epetra_Vector> densityDotn,
      Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);

    //@}

  protected:

    /// calculate interaction with neighboring particles
    virtual void CalcNeighboringParticleMeshFreeInteraction(
      const int lidNodeCol_i,
      const std::list<DRT::Node*> neighboring_particles,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<Epetra_Vector> densityDotn,
      Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);

    /// calculate interaction with neighboring walls
    virtual void CalcNeighboringWallMeshFreeInteraction(
      const int lidNodeCol_i,
      const std::set<DRT::Element*>& neighboring_walls,
      const Teuchos::RCP<DRT::Discretization>& walldiscret,
      const Teuchos::RCP<const Epetra_Vector>& walldisn,
      const Teuchos::RCP<const Epetra_Vector>& wallveln,
      const Teuchos::RCP<Epetra_Vector>& accn,
      Teuchos::RCP<Epetra_Vector> densityDotn);

    virtual void CalcNeighboringHeatSourcesContact(
      const int lidNodeCol_i,
      const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources,
      const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn);

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

    /// temporary storage of wallMeshFree data
    WallMeshFreeData wallMeshFreeData_;

    //@}

    //! @name wall interaction parameters
    //@{

    /// laplacian velocity, diffusion coefficient
    double diffusionCoeff_;

    /// gradient of the divergence of the velocity, convection coefficient
    double convectionCoeff_;

    //@}

  };  // class ParticleCollisionHandler

  //! @name Container structs for faster access
  //@{

  /// struct to store particleMeshFree data for faster access
  struct ParticleMeshFreeData {
    LINALG::Matrix<3,1> dis;
    LINALG::Matrix<3,1> vel;
/// --- second round -> unactivated ----------------------------------- ///
//    LINALG::Matrix<3,1> gradT_over_rho;
///--------------------------------------------------------------------///
    double radius;
    double density;
    double specEnthalpy;
    double mass;
    double pressure;
    double temperature;
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
