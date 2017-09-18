/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_interaction.H

\brief Particle-MeshFree interaction handling

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

/*
 References:
 Antoci2007: Numerical simulation of fluid–structure interaction by SPH, doi: 10.1016/j.compstruc.2007.01.002.
 Monaghan2005: Smoothed particle hydrodynamics, doi:10.1088/0034-4885/68/8/R01.
 Espanol2003: Smoothed dissipative particle dynamics, DOI: 10.1103/PhysRevE.67.026705.
 Akinci2013: http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf
 Adami2012: A generalized wall boundary condition for smoothed particle hydrodynamics, http://dx.doi.org/10.1016/j.jcp.2012.05.005
 Adami2013: A transport-velocity formulation for smoothed particle hydrodynamics, http://dx.doi.org/10.1016/j.jcp.2013.01.043
*/

#ifndef SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_
#define SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_
/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle.H"
#include "particle_node.H"
#include "particle_interData.H"
#include "particle_particle.H"
#include "particleMeshFree_weightFunction.H"
#include "../drt_lib/drt_element.H"
#include <boost/unordered_map.hpp>
#include "particle_enums.H"
#include "../linalg/linalg_sparsematrix.H"

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

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

namespace LINALG
{
  class SparseMatrix;
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
      const double initDensity,
      const double restDensity,
      const double refdensfac);

    //! Destructor
    virtual ~ParticleMeshFreeInteractionHandler(){ ; }

    //@}

    //! @name Actions
    //@{

    //! set up the object for the interaction computations
    virtual void Init(
        const int step,
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn);

    virtual void Init(
        const int step,
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> temperature);

    virtual void Init(
        const int step,
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> mass,
        Teuchos::RCP<const Epetra_Vector> specEnthalpyn,
        Teuchos::RCP<const Epetra_Vector> temperature,
        Teuchos::RCP<const Epetra_Vector> densityn,
        Teuchos::RCP<const Epetra_Vector> pressure);

    //! set up the ColParticles
    virtual void InitColParticles();

    //! accelerations, modified pressures and velocities for boundary particles
    virtual void InitBoundaryData(
        Teuchos::RCP<const Epetra_Vector> accn,
        const LINALG::Matrix<3,1>& g,
        double &bpintergy);

    //! add to the neighbour maps the neighbours of the present time step
    virtual void AddNewNeighbours(const int step);

    //! update weights and distances in all the neighbours (all maps, all the particles, all the walls)
    //! to speed-up the process the interactions with "step" are skipped. -1 means that nothing is skipped
    virtual void UpdateWeights(const int step = -1);

    //! update neighbour weights

    //! set state vector in ColParticles
    virtual void SetStateVector(
        Teuchos::RCP<const Epetra_Vector> colorFieldGradientn,
        const PARTICLE::StateVectorType svt);

    //! set state vector in ColParticles - for multiVector
    //virtual void SetStateVector(Teuchos::RCP<const Epetra_MultiVector> stateVector);

    //! clear content, keep memory
    virtual void Clear();

    //! clear content, keep memory and keep interactions up to step-memory
    virtual void Clear(const int step, const int memory);

    //! Apply coordinate shift in case of periodic boundary conditions
    virtual void ShiftPeriodicBoundaryPair(
        const LINALG::Matrix<3,1>& dis_i,
        LINALG::Matrix<3,1>& dis_j,
        const double& radius_1,
        const double& radius_2);

    //@}

    //! @name interactions
    //@{

    /// compute the density part of the particle vs particle interactions
    virtual void Inter_pvp_densityDot(
        const Teuchos::RCP<Epetra_Vector> densityDotn,
        const double extMulti = 1.0);

    /// compute the density part of the boundary particle vs particle interactions
    virtual void Inter_bpvp_densityDot(
        const Teuchos::RCP<Epetra_Vector> densityDotn,
        const double extMulti = 1.0);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_pvp_acc_var1(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> p0_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double extMulti = 1.0);

    /// compute the momentum part of the boundary particle vs particle interactions
    virtual void Inter_bpvp_acc_var1(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> p0_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double extMulti = 1.0);

    /// compute the momentum part of the particle vs particle interactions: alternative variant according to Adami2013
    virtual void Inter_pvp_acc_var2(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> p0_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double extMulti = 1.0);

    /// compute the momentum part of the boundary particle vs particle interactions: alternative variant according to Adami2013
    virtual void Inter_bpvp_acc_var2(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> p0_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double extMulti = 1.0);

    ///Computes Gradient of divergence operator for vector v - pvp
    virtual void Inter_pvp_Laplace_x(
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij,
        const LINALG::Matrix<3,1>& x_ij,
        LINALG::Matrix<3,1>& Laplace_x_ij);

    /// Computes Laplace operator for vector v - pvp
    virtual void Inter_pvp_GradDiv_x(
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij,
        const LINALG::Matrix<3,1>& x_ij,
        LINALG::Matrix<3,1>& GradDiv_x_ij);

    /// compute the heat of the particle vs particle interactions
    virtual void Inter_pvp_specEnthalpyDot(
        const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
        const double extMulti = 1.0);

    /// compute the colorFieldGradient of the particle vs particle interactions
    virtual void Inter_pvp_colorFieldGradient(
        const Teuchos::RCP<Epetra_Vector> colorFieldGradientn,
        const double extMulti = 1.0);

    /// compute the accn from the color field gradient of the particle vs particle interactions
    virtual void Inter_pvp_surfaceTensionCFG(
        const Teuchos::RCP<Epetra_Vector> accn,
        const double extMulti = 1.0);

    /// compute the gradient of the momentum part of the particle vs particle interactions
    /// if you change this function you need to touch the pressure definition too!
    virtual void Inter_pvp_gradAccP(
        const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
        const double &restDensity,
        const double extMulti = 1.0);

    virtual void Inter_pvp_gradAccPapproxOnlyHess(
        const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
        const double &restDensity,
        const double extMulti = 1.0);

    /// compute the density part of the particle vs wall interactions
    virtual void Inter_pvw_densityDot(
        const Teuchos::RCP<Epetra_Vector> densityDotn,
        const double extMulti = 1.0);

    /// compute the momentum part of the particle vs wall interactions
    virtual void Inter_pvw_acc(
        const Teuchos::RCP<Epetra_Vector> accn,
        const double extMulti = 1.0);

    /// compute the heat of the particle vs heat sources interactions
    virtual void Inter_pvhs_specEnthalpyDot(
        const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
        const double extMulti = 1.0);

    virtual void Inter_pvw_gradAccP(
        const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
        const double &restDensity,
        const double extMulti = 1.0);

    virtual void Inter_pvw_gradAccPapproxOnlyHess(
        const Teuchos::RCP<LINALG::SparseMatrix> gradAccP,
        const double &restDensity,
        const double extMulti = 1.0);

    //@}

    //! @name mesh free fields
    //@{

    /// compute density - mesh free sense
    virtual void MF_mW(
        const Teuchos::RCP<Epetra_Vector> mW,
        const Teuchos::RCP<Epetra_Vector> unity_vec = Teuchos::null,
        const double extMulti = 1.0);

    /// compute density - mesh free sense - wall contributions
    virtual void MF_mW_Wall(
        const Teuchos::RCP<Epetra_Vector> mW,
        const double extMulti = 1.0);

    /// compute density - mesh free sense - boundary particle contributions
    virtual void MF_mW_Boundary(
        const Teuchos::RCP<Epetra_Vector> mW,
        const Teuchos::RCP<Epetra_Vector> unity_vec,
        const double extMulti = 1.0);

    /// re-initialize density in case of free-surface flow - mesh free sense
    virtual void MF_ReInitDensity(
        const Teuchos::RCP<Epetra_Vector> density,
        const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType);

    /// compute pressure difference (for debugging)
    virtual void MF_DeltaP(
        const Teuchos::RCP<Epetra_Vector> PressureDiff);

    /// initialize free surface particles
    virtual void InitFreeSurfaceParticles();

    /// set pressure of free surface particles
    virtual void SetPressureOfFreeSurfaceParticles(
        const Teuchos::RCP<Epetra_Vector> densityn,
        const Teuchos::RCP<Epetra_Vector> pressuren);

    /// determine accelerations of free-surface particle
    virtual void Inter_fsvp_acc(
        const Teuchos::RCP<Epetra_Vector> acc);

    /// compute m * gradW - mesh free sense
    virtual void MF_mGradW(
        const Teuchos::RCP<Epetra_Vector> mGradW,
        const bool withWalls = true,
        const double extMulti = 1.0);

    /// compute m * hessW - mesh free sense
    virtual void MF_mHessW(
        const Teuchos::RCP<Epetra_MultiVector> mHessW,
        const bool withWalls = true,
        const double extMulti = 1.0);

    /// Compute linearizaton of AccP via FAD (for debugging) - mesh free style
    void Inter_pvp_gradAccP_FAD(
        const Teuchos::RCP<Epetra_Vector> res,
        const Teuchos::RCP<LINALG::SparseMatrix> stiff);

    /// Compute consistent linearizaton of AccP - mesh free style
    void Inter_pvp_gradAccP_consistent(
        const Teuchos::RCP<Epetra_Vector> res,
        const Teuchos::RCP<LINALG::SparseMatrix> stiff);

    //@}

    //! @name mesh free debug tools
    //@{

    /// wrapper for the template
    void PrintNeighbours(const PARTICLE::NeighbourType& nt);

    /// Get min / max pressure etc.
    void GetExtremeValues(
        double& min_pvp_dist,
        double& min_pvw_dist,
        double& min_pressure,
        double& max_pressure)
    {
      min_pvp_dist=min_pvp_dist_;
      min_pvw_dist=min_pvw_dist_;
      min_pressure=min_pressure_;
      max_pressure=max_pressure_;
    }

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

    //! add to the neighbour maps the neighbours of the present time step - particles (overlapping vector)
    virtual void AddNewNeighbours_op(
        const ParticleMF& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p,
        const int step);

    //! add to the neighbour maps the neighbours of the present time step - boundary particles
    virtual void AddNewNeighbours_bp(
        const ParticleMF& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p,
        const int step);

    //! add to the neighbour maps the neighbours of the present time step - free-surface particles
    virtual void AddNewNeighbours_fsp(
        const ParticleMF& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p,
        const bool mark_fs_indirect=false,
        const int step=1);

    //! add to the neighbour maps the neighbours of the present time step - walls
    virtual void AddNewNeighbours_w(
        const ParticleMF& particle_i,
        const boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
        const int step);

    //! add to the neighbour maps the neighbours of the present time step - heat sources
    virtual void AddNewNeighbours_hs(
        const ParticleMF& particle_i,
        const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > >& neighboring_heatSources);

    //! update the neighbours_p_ map
    virtual void UpdateWeights_p(const int step = -1);

    //! update the neighbours_w_ map
    virtual void UpdateWeights_w(const int step = -1);

    //! update weights and distances of boundary particle neighbours
    virtual void UpdateWeights_bp(const int step = -1);

     //! update the neighbours_op_ map
     virtual void UpdateWeights_op(const int step = -1);

     //! build periodic boundary conditions
     virtual void BuildPeriodicBC();

    /// template to print neighbours
    template <typename obj>
    void PrintNeighbours(const std::vector< boost::unordered_map<int, obj> >& neighbours)
    {
      std::cout << "\n\n interactions\n\n";
      //bool trg_interactions = false;
      for (unsigned int lidNodeRow_i = 0; lidNodeRow_i != neighbours.size(); ++lidNodeRow_i)
      {
        const int gid_i = discret_->NodeRowMap()->GID(lidNodeRow_i);

        std::cout << "gid " << gid_i << " | ";

        typename boost::unordered_map<int, obj >::const_iterator jj;
        for (jj = neighbours[lidNodeRow_i].begin(); jj != neighbours[lidNodeRow_i].end(); ++jj)
        {
          std::cout << " " << jj->first;
          //trg_interactions = true;
        }
        std::cout << std::endl;
      }
    }

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

    /// temporary storage for FAD version of particleMeshFree data
    std::vector<Teuchos::RCP<ParticleFAD> > colFADParticles_;

    /// neighbours - particles
    std::vector<boost::unordered_map<int, InterDataPvP> > neighbours_p_;

    /// neighbours - particles (overlapping storage)
    std::vector<boost::unordered_map<int, InterDataPvP> > overlappingneighbours_p_;

    /// neighbours - heat sources
    std::vector<boost::unordered_map<int, Teuchos::RCP<HeatSource> > > neighbours_hs_;

    /// neighbours - walls
    std::vector<boost::unordered_map<int, InterDataPvW> > neighbours_w_;

    /// weight function
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    /// temporary storage of wallMeshFree data
    WallMeshFreeData wallMeshFreeData_;

    /// neighbours - boundary particles
    std::map<int,std::map<int, InterDataPvP> > neighbours_bp_;

    /// neighbours - free surface particles
    std::map<int,std::map<int, InterDataPvP> > neighbours_fsp_;

    /// standard map that contains pointers to all boundary particles of column map
    std::map<int,ParticleMF* > boundaryparticles_;

    /// standard map that contains pointers to all free-surface particles of column map
    std::map<int,ParticleMF* > freesurfaceparticles_;

    /// number of space dimensions required for weight function definition
    const INPAR::PARTICLE::WeightFunctionDim WF_DIM_;

    //! rest density for meshFree-pressure dynamics
    const double initDensity_;

    //! rest density for meshFree-pressure dynamics defining the equation of state: pressure=c_0^2(density-refdensfac_*restDensity_)  with the speed of sound c_0^2=BULK_MODULUS/restDensity_
    const double restDensity_;

    /// factor defining the equation of state: pressure=c_0^2(density-refdensfac_*restDensity_)  with the speed of sound c_0^2=BULK_MODULUS/restDensity_;
    const double refdensfac_;

    /// laplacian velocity, diffusion coefficient
    double diffusionCoeff_;

    /// gradient of the divergence of the velocity, convection coefficient
    double convectionCoeff_;

    /// artificial viscosity according to Adami et al. 2012, Eq. (11)
    double artificialViscosity_;

    /// cohesion surface tension coefficient with void (or air)
    double surfaceVoidTension_;

    /// cohesion surface tension coefficient with the wall
    double surfaceWallTension_;

    // security trigger in case the color field gradient was not set properly
    bool trg_updatedColorFieldGradient_;

    // length of periodic box if appicable (=-1.0 if not)
    double periodic_length_;

    // minimal distance (over all particles) between two particles
    double min_pvp_dist_;

    // minimal distance (over all particles) between particles and walls
    double min_pvw_dist_;

    // minimal pressure (over all particles)
    double min_pressure_;

    // maximal pressure (over all particles)
    double max_pressure_;

    //@}

  };  // class ParticleCollisionHandler

}  // namespace PARTICLE
/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_INTERACTION_H_ */
