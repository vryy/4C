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
#include "particle_node.H"
#include "particle_interData.H"
#include "particle_particle.H"
#include "particleMeshFree_weightFunction.H"
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

    //! set up the ColParticles
    virtual void InitColParticles();

    //! accelerations, modified pressures and velocities for boundary particles
    virtual void InitBoundaryData(
        Teuchos::RCP<const Epetra_Vector> accn,
        const LINALG::Matrix<3,1>& g,
        double &bpintergy);

    //! add to the neighbour maps the neighbours of the present time step
    virtual void AddNewNeighbours(const int step);

    //! set state vector in ColParticles
    virtual void SetStateVector(
        Teuchos::RCP<const Epetra_Vector> stateVector,
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
    virtual void Inter_pvp_acc(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> trvl_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double time = -1.0);

    /// compute the momentum part of the boundary particle vs particle interactions
    virtual void Inter_bpvp_acc(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> trvl_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double time = -1.0);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_fspvp_acc1(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> curvature,
        double time);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_fspvp_acc1b(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> curvature,
        double time);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_fspvp_acc2(
        const Teuchos::RCP<Epetra_Vector> accn,
        double time);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_fspvp_acc2b(
        const Teuchos::RCP<Epetra_Vector> accn,
        double time);

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

    /// compute the heat of the particle vs heat sources interactions
    virtual void Inter_pvhs_specEnthalpyDot(
        const Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
        const double extMulti = 1.0);

    //@}

    //! @name mesh free fields
    //@{

    /// compute density - mesh free sense
    virtual void MF_mW(
        const Teuchos::RCP<Epetra_Vector> mW,
        const Teuchos::RCP<Epetra_Vector> unity_vec = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> unity_vec_grad = Teuchos::null,
        const double extMulti = 1.0);

    /// compute density - mesh free sense - boundary particle contributions
    virtual void MF_mW_Boundary(
        const Teuchos::RCP<Epetra_Vector> mW,
        const Teuchos::RCP<Epetra_Vector> unity_vec,
        const Teuchos::RCP<Epetra_Vector> unity_vec_grad,
        const double extMulti = 1.0);

    void MF_SmoothedCFG(
        const Teuchos::RCP<Epetra_Vector> unity_vec_grad);

    /// re-initialize density in case of free-surface flow - mesh free sense
    virtual void MF_ReInitDensity(
        const Teuchos::RCP<Epetra_Vector> density,
        const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType);

    /// initialize free surface particles
    virtual void InitFreeSurfaceParticles(
        const Teuchos::RCP<Epetra_Vector> fspType);

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

    //! @name AddNewNeighbours - these functions are collected in AddNewNeighbours for faster search
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

    //! add to the neighbour maps the neighbours of the present time step - heat sources
    virtual void AddNewNeighbours_hs(
        const ParticleMF& particle_i,
        const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > >& neighboring_heatSources);

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

    //! @name interaction methods
    //@{

    //! construct the specific coefficients
    virtual void Inter_generalCoeff(
        double& generalCoeff_ij,
        double& generalCoeff_ji,
        const double& density_i,
        const double& mass_i,
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate pressure
    virtual void Inter_pressure(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ji,
        const double& generalCoeff_ij,
        const double& generalCoeff_ji,
        const double& density_i,
        const double& pressure_i,
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate background pressure
    virtual void Inter_backgroundPressure(
        LINALG::Matrix<3,1>* trvl_accn_ij,
        LINALG::Matrix<3,1>* trvl_accn_ji,
        const double& generalCoeff_ij,
        const double& generalCoeff_ji,
        const double& density_i,
        const double& mass_i,
        const double& radius_i,
        const double& pressure_i,
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate xsph contribution
    virtual void Inter_xsph(
        LINALG::Matrix<3,1>* trvl_accn_ij,
        LINALG::Matrix<3,1>* trvl_accn_ji,
        const double& density_i,
        const double& mass_i,
        const double& radius_i,
        const ParticleMF& particle_j,
        const LINALG::Matrix<3,1>& vRel_ij,
        const InterDataPvP& interData_ij
        );

    //! evaluate additional term divA due to transport velocity
    virtual void Inter_transportVelocity_divA(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ji,
        LINALG::Matrix<3,1>* accn_ij_A,
        LINALG::Matrix<3,1>* accn_ji_A,
        const double& generalCoeff_ij,
        const double& generalCoeff_ji,
        const double& density_i,
        const LINALG::Matrix<3,1>& vel_i,
        const LINALG::Matrix<3,1>& velConv_i,
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate laminar viscosity
    virtual void Inter_laminarViscosity(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ji,
        const double& generalCoeff_ij,
        const double& generalCoeff_ji,
        const double& density_i,
        const ParticleMF& particle_j,
        const LINALG::Matrix<3,1>& vRel_ij,
        const InterDataPvP& interData_ij
        );

    //! evaluate artificial viscosity
    virtual void Inter_artificialViscosity(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ji,
        const double& density_i,
        const double& mass_i,
        const double& radius_i,
        const ParticleMF& particle_j,
        const LINALG::Matrix<3,1>& vRel_ij,
        const InterDataPvP& interData_ij
        );

    //! evaluate surface tension
    virtual void Inter_surfaceTension(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ji,
        const double& density_i,
        const double& mass_i,
        const double& radius_i,
        const ParticleMF& particle_j,
        const InterDataPvP& interData_ij,
        const double& time
        );

    //! evaluate surface tension interaction potential
    virtual void SurfTensionInterPot(
        const double& radius,
        const double& rRelNorm2,
        double& lambda,
        double& potential
    );

    //! return surface tension time fac
    virtual double SurfTensionTimeFac(
        const double& time
    );

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

    /// weight function
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    /// particle interaction variant
    const bool interactionVariant2_;

    /// solve thermal problem
    const bool solve_thermal_prob_;

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    /// free surface type
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType_;

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

    /// background pressure employed for modified particle convection velocity
    const double background_pressure_;

    /// apply particle convection velocity that differs from momentum velocity
    const bool transport_velocity_;

    /// do not apply velocity difference tensor in case of transport velocity formulation
    const bool no_veldiff_term_;

    /// apply additional damping force in order to achieve static equilibrium solution
    const double damping_factor_;

    /// damping factor of XSPH scheme
    const double xsph_dampfac_;

    /// stiffness factor of XSPH scheme
    const double xsph_stiffac_;

    /// laplacian velocity, diffusion coefficient
    double diffusionCoeff_;

    /// gradient of the divergence of the velocity, convection coefficient
    double convectionCoeff_;

    /// artificial viscosity according to Adami et al. 2012, Eq. (11)
    double artificialViscosity_;

    /// surface tension between fluid and void phase
    double surfaceTension_;

    /// type of surface tension model
    const INPAR::PARTICLE::SurfaceTensionType surfaceTensionType_;

    /// static contact angle between fluid and solid phase in degree measure
    double staticContactAngle_;

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
