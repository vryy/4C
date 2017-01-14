/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.H

\brief Particle-MeshFree interaction handling

\level 3

\maintainer Alessandro Cattabiani
*/

/*
 References:
 Antoci2007: Numerical simulation of fluid–structure interaction by SPH, doi: 10.1016/j.compstruc.2007.01.002.
 Monaghan2005: Smoothed particle hydrodynamics, doi:10.1088/0034-4885/68/8/R01.
 Espanol2003: Smoothed dissipative particle dynamics, DOI: 10.1103/PhysRevE.67.026705.
 Akinci2013: http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf
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
  struct Interaction_Particles;
  struct Interaction_Walls;

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
      const Teuchos::ParameterList& particledynparams,
      const double restDensity);

    //! Destructor
    virtual ~ParticleMeshFreeInteractionHandler(){ ; }

    //@}

    //! @name Actions
    //@{

    //! set up the object for the interaction computations
    virtual void Init(
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> densityapproxn,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> temperature,
        Teuchos::RCP<const Epetra_Vector> pressure);

    /// clear content, keep memory
    virtual void Clear();

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual void EvaluateParticleMeshFreeInteractions(
      const Teuchos::RCP<Epetra_Vector> accn,
      const Teuchos::RCP<Epetra_Vector> densityDotn,
      const Teuchos::RCP<Epetra_Vector> densityapproxn,
      const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);

    //@}

  protected:

    //! set state vectors (col type)
    virtual void SetStateVectors(
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> densityapproxn,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> temperature,
        Teuchos::RCP<const Epetra_Vector> pressure);

    //! set state vectors (col type) - second round
    virtual void SetStateVectors_SecondRound(
        Teuchos::RCP<const Epetra_Vector> colorFieldGradientn);

    /// determine the neighbors
    virtual void SetNeighbours();

    // determine the neighbors - particles
    void SetNeighbours_Particles(
        const ParticleMeshFreeData& particle_i,
        const std::list<DRT::Node*>& neighboring_particles);

    // determine the neighbors - walls
    void SetNeighbours_Walls(
        const ParticleMeshFreeData& particle_i,
        const std::set<DRT::Element*>& neighboring_walls,
        const Teuchos::RCP<DRT::Discretization>& walldiscret,
        const Teuchos::RCP<const Epetra_Vector>& walldisn,
        const Teuchos::RCP<const Epetra_Vector>& wallveln);

    // determine the neighbors - heat sources
    void SetNeighbours_HeatSources(
        const ParticleMeshFreeData& particle_i,
        const std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less>& neighboring_heatSources);

    // print interactions
    void PrintNeighbouringParticles();

    /// calculate interaction with neighboring particles
    virtual void CalcNeighboringParticleMeshFreeInteraction(
      const ParticleMeshFreeData& particle_i,
      const std::list<Interaction_Particles> neighboringParticles,
      const Teuchos::RCP<Epetra_Vector> accn,
      const Teuchos::RCP<Epetra_Vector> densityDotn,
      Teuchos::RCP<Epetra_Vector> densityapproxn,
      const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
      const Teuchos::RCP<Epetra_Vector> colorFieldGradientn);

    /// calculate interaction with neighboring walls
    virtual void CalcNeighboringWallMeshFreeInteraction(
      const ParticleMeshFreeData& particle_i,
      const std::list<Interaction_Walls>& neighbouringWalls,
      const Teuchos::RCP<Epetra_Vector> accn,
      const Teuchos::RCP<Epetra_Vector> densityDotn);

    /// calculate interaction with heat sources
    virtual void CalcNeighboringHeatSourcesContact(
      const ParticleMeshFreeData& particle_i,
      const std::list<Teuchos::RCP<HeatSource> > neighboringHeatSources,
      const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn);

    /// calculate interaction with neighboring particles - second round
    virtual void CalcNeighboringParticleMeshFreeInteraction_SecondRound(
      const ParticleMeshFreeData& particle_i,
      const std::list<Interaction_Particles> neighbouringParticles,
      const Teuchos::RCP<Epetra_Vector> accn);

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

    /// neighbours - particles
    std::vector<std::list<Interaction_Particles> > neighbouringParticles_;

    /// neighbours - heat sources
    std::vector<std::list<Teuchos::RCP<HeatSource> > > neighbouringHeatSources_;

    /// neighbours - heat sources
    std::vector<std::list<Interaction_Walls> > neighbouringWalls_;

    /// weight function derivative pointer
    double (*weightDerivative_)(const double &, const double &);

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    /// temporary storage of wallMeshFree data
    WallMeshFreeData wallMeshFreeData_;

    //@}

    //! @name wall interaction parameters
    //@{

    /// rest density
    const double restDensity_;

    /// laplacian velocity, diffusion coefficient
    double diffusionCoeff_;

    /// gradient of the divergence of the velocity, convection coefficient
    double convectionCoeff_;

    /// cohesion surface tension coefficient with void (or air)
    double surfaceVoidTension_;

    /// cohesion surface tension coefficient with the wall
    double surfaceWallTension_;

    // trigger for the second round
    bool trg_secondRound_;

    //@}

  };  // class ParticleCollisionHandler

  //! @name Container structs for faster access
  //@{

  /// struct to store particleMeshFree data for faster access
  struct ParticleMeshFreeData {
    LINALG::Matrix<3,1> dis;
    LINALG::Matrix<3,1> vel;
    LINALG::Matrix<3,1> colorFieldGradient; //< Between E2 and E3 in http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf

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

  /// struct to store some basic interaction data pvp
  struct Interaction_Particles {
    int lidCol_j;
    LINALG::Matrix<3,1> rRelVersor_ij;
    double weightDerivative_ij;
    double weightDerivative_ji;
    double rRelNorm2;
  };

  /// struct to store some basic interaction data pvw
  struct Interaction_Walls {
    int elemId;
    LINALG::Matrix<3,1> rRelVersor;
    double weightDerivative;
    double rRelNorm2;
  };

  //@}

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
