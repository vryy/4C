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
#include "particle_interData.H"
#include "particle_particle.H"
#include "particleMeshFree_weightFunction.H"
#include "../drt_lib/drt_element.H"
#include <boost/unordered_map.hpp>

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

  /// struct to store wallMeshFree data for faster access
  struct WallMeshFreeData {
    double density_;
    double mass_;
    double pressure_;
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
  enum StateVectorType
  {
    Dis,
    Vel,
    ColorFieldGradient,
    Radius,
    Density,
    DensityDot,
    SpecEnthalpy,
    Mass,
    Pressure,
    Temperature,
    Alpha
  };

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
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> pressure,
        Teuchos::RCP<const Epetra_Vector> temperature,
        Teuchos::RCP<const Epetra_Vector> densityapproxn,
        const int step);

    //! set up the ColParticles
    virtual void InitColParticles(
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> pressure,
        Teuchos::RCP<const Epetra_Vector> temperature,
        Teuchos::RCP<const Epetra_Vector> densityapproxn);

    //! add to the neighbour maps the neighbours of the present time step
    virtual void AddNewNeighbours(const int step);

    //! update weights and distances in all the neighbours (all maps, all the particles, all the walls)
    //! to speed-up the process the interactions with "step" are skipped. -1 means that nothing is skipped
    virtual void UpdateWeights(const int step = -1) { UpdateWeights_p(step); UpdateWeights_w(step); }

    //! update neighbour weights

    //! set state vector in ColParticles - color field gradient
    virtual void SetStateVector(
        Teuchos::RCP<const Epetra_Vector> colorFieldGradientn,
        const PARTICLE::StateVectorType svt);

    //! clear content, keep memory
    virtual void Clear();

    //! clear content, keep memory and keep interactions up to step-memory
    virtual void Clear(const int step, const int memory);

    //@}

    //! @name interactions
    //@{

    /// compute the density part of the particle vs particle interactions
    virtual void Inter_pvp_densityDot(const Teuchos::RCP<Epetra_Vector> densityDotn);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_pvp_acc(const Teuchos::RCP<Epetra_Vector> accn, const bool trg_pressure = true);

    /// compute the heat of the particle vs particle interactions
    virtual void Inter_pvp_specEnthalpyDot(const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);

    /// compute the colorFieldGradient of the particle vs particle interactions
    virtual void Inter_pvp_colorFieldGradient(const Teuchos::RCP<Epetra_Vector> colorFieldGradientn);

    /// compute the accn from the color field gradient of the particle vs particle interactions
    virtual void Inter_pvp_surfaceTensionCFG(const Teuchos::RCP<Epetra_Vector> accn);

    /// compute the density part of the particle vs wall interactions
    virtual void Inter_pvw_densityDot(const Teuchos::RCP<Epetra_Vector> densityDotn);

    /// compute the momentum part of the particle vs wall interactions
    virtual void Inter_pvw_acc(const Teuchos::RCP<Epetra_Vector> accn);

    /// compute the heat of the particle vs heat sources interactions
    virtual void Inter_pvhs_specEnthalpyDot(const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn);

    /// compute the divergence-free pressure forces adjusting the stiffness parameters
    /// Steps 6 and 7 (only the third term) of Algorithm 2 in http://dx.doi.org/10.1145/2786784.2786796
    virtual void Inter_pvp_divFreePressureAcc(
        const Teuchos::RCP<Epetra_Vector> divFreePressureAcc,
        const double dt);

    /// compute the constant-density pressure forces adjusting the stiffness parameters
    /// Steps 6 and 7 (only the third term) of Algorithm 3 in http://dx.doi.org/10.1145/2786784.2786796
    virtual void Inter_pvp_constDensityPressureAcc(
        const Teuchos::RCP<Epetra_Vector> divFreePressureAcc,
        const double dt);

    //@}

    //! @name mesh free fields
    //@{

    /// compute density and update colParticles - mesh free sense
    virtual void MF_density(const Teuchos::RCP<Epetra_Vector> densityn);

    /// compute alpha factors and update colParticles - mesh free sense
    // http://dx.doi.org/10.1145/2786784.2786796
    // Equation 8
    virtual void MF_alpha();

    //@}

    //! @name print methods
    //@{

    // print interactions
    void Print_neighbours_p();

    //@}

  protected:

    //! @name AddNewNeighbours - these functions are collected in AddNewNeighbours for faster search. This is the most modular
    //! way to compute neighbours that I can think off (without hindering performances) If you have a better idea I am open to suggestions
    //@{

    //! add to the neighbour maps the neighbours of the present time step - particles
    virtual void AddNewNeighbours_p(
        const ParticleMF& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p,
        const int step);

    //! add to the neighbour maps the neighbours of the present time step - walls
    virtual void AddNewNeighbours_w(
        const ParticleMF& particle_i,
        const boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
        const Teuchos::RCP<DRT::Discretization>& walldiscret,
        const Teuchos::RCP<const Epetra_Vector>& walldisn,
        const Teuchos::RCP<const Epetra_Vector>& wallveln,
        const int step);

    //! add to the neighbour maps the neighbours of the present time step - heat sources
    virtual void AddNewNeighbours_hs(
        const ParticleMF& particle_i,
        const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > >& neighboring_heatSources);

    virtual void UpdateWeights_p(const int step = -1);

    virtual void UpdateWeights_w(const int step = -1);

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
    std::vector<ParticleMF> colParticles_;

    /// neighbours - particles
    std::vector<boost::unordered_map<int, InterDataPvP> > neighbours_p_;

    /// neighbours - heat sources
    std::vector<boost::unordered_map<int, Teuchos::RCP<HeatSource> > > neighbours_hs_;

    /// neighbours - heat sources
    std::vector<boost::unordered_map<int, InterDataPvW> > neighbours_w_;

    /// weight function
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    /// temporary storage of wallMeshFree data
    WallMeshFreeData wallMeshFreeData_;

    //@}

    //! @name wall interaction parameters
    //@{

    /// rest density
    const double restDensity_;

    /// lower cap of the alpha parameter
    const double alphaMin_;

    /// laplacian velocity, diffusion coefficient
    double diffusionCoeff_;

    /// gradient of the divergence of the velocity, convection coefficient
    double convectionCoeff_;

    /// cohesion surface tension coefficient with void (or air)
    double surfaceVoidTension_;

    /// cohesion surface tension coefficient with the wall
    double surfaceWallTension_;

    // security trigger in case the color field gradient was not set properly
    bool trg_updatedColorFieldGradient_;

    //@}

  };  // class ParticleCollisionHandler

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
