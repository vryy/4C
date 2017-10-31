/*----------------------------------------------------------------------*/
/*!
\file particle_sph_interaction.H

\brief Particle SPH interaction handling

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

#ifndef SRC_DRT_PARTICLE_PARTICLE_SPH_INTERACTION_H_
#define SRC_DRT_PARTICLE_PARTICLE_SPH_INTERACTION_H_
/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_discret.H"
#include "particle_interData.H"
#include "particle_particle.H"
#include <boost/unordered_map.hpp>
#include "particle_sph_weightFunction.H"

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

  /// struct to store wallSPH data for faster access
  struct WallSPHData {
    double density_;
    double mass_;
    double pressure_;
  };
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
   * \brief Interaction handler for particle SPH interactions
   *
   * \author meier
   * \date 10/17
   */
  class ParticleSPHInteractionHandler
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    ParticleSPHInteractionHandler(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
      const Teuchos::ParameterList& particledynparams,
      bool norender=false);

    //! Destructor
    virtual ~ParticleSPHInteractionHandler(){ ; }

    //@}

    //! @name Actions
    //@{

    //! set up the object for the interaction computations
    virtual void Init(
        Teuchos::RCP<const Epetra_Vector> disn,
        Teuchos::RCP<const Epetra_Vector> veln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<const Epetra_Vector> mass);

    //! set up the ColParticles
    virtual void InitColParticles();

    //! set density & mass at the very beginning
    void InitDensityAndMass(
        const double particle_volume,
        const Teuchos::RCP<Epetra_Vector> density,
        const Teuchos::RCP<Epetra_Vector> mass);

    //! accelerations, modified pressures and velocities for boundary particles
    virtual void InitBoundaryData(
        Teuchos::RCP<const Epetra_Vector> accn,
        const LINALG::Matrix<3,1>& g);

    //! add to the neighbour maps the current neighbours
    virtual void AddNewNeighbours();

    //! set state vector in ColParticles
    virtual void SetStateVector(
        Teuchos::RCP<const Epetra_Vector> stateVector,
        const PARTICLE::StateVectorType svt);

    //! clear content, keep memory
    virtual void Clear();

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
        const Teuchos::RCP<Epetra_Vector> densityDotn);

    /// compute the momentum part of the particle vs particle interactions
    virtual void Inter_pvp_acc(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> trvl_acc = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::null,
        const double time = -1.0);

    /// compute surface tension interaction according to continuum model by Adami et al. 2010
    virtual void Inter_fspvp_Adami_1(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> curvature,
        double time);

    /// compute surface tension interaction according to continuum model by Adami et al. 2010 (modified version)
    virtual void Inter_fspvp_Adami_2(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> curvature,
        double time);

    /// compute surface tension interaction according to continuum model by Adami et al. 2010 (modified version)
    virtual void Inter_fspvp_Adami_3(
        const Teuchos::RCP<Epetra_Vector> accn,
        const Teuchos::RCP<Epetra_Vector> curvature,
        double time);

    /// compute surface tension interaction according to continuum model by Hu et al. 2006
    virtual void Inter_fspvp_Hu(
        const Teuchos::RCP<Epetra_Vector> accn,
        double time);

    /// compute surface tension interaction according to continuum model by Hu et al. 2006 (modified version)
    virtual void Inter_fspvp_Hu_b(
        const Teuchos::RCP<Epetra_Vector> accn,
        double time);

    ///Computes Gradient of divergence operator for vector v - pvp
    virtual void Inter_pvp_Laplace_x(
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij,
        const LINALG::Matrix<3,1>& x_ij,
        LINALG::Matrix<3,1>& Laplace_x_ij);

    /// Computes Laplace operator for vector v - pvp
    virtual void Inter_pvp_GradDiv_x(
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij,
        const LINALG::Matrix<3,1>& x_ij,
        LINALG::Matrix<3,1>& GradDiv_x_ij);

    //! Compute pressure for the entire vector
    virtual void Density2Pressure(
        const Teuchos::RCP<const Epetra_Vector> density,
        Teuchos::RCP<Epetra_Vector> &pressure);

    //! Determine internal energy
    virtual void DetermineIntEnergy(double &intenergy);

    //@}

    //! @name mesh free fields
    //@{

    /// compute \sum m * W (usually the density)
    virtual void MF_mW(
        const Teuchos::RCP<Epetra_Vector> mW,
        const Teuchos::RCP<Epetra_Vector> unity_vec = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> unity_vec_grad = Teuchos::null,
        const Teuchos::RCP<Epetra_Vector> phaseColor = Teuchos::null);

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

    /// Determine coefficients required for viscous forces
    void ViscousCoefficients(const double viscosity, const double bulkViscosity, double& convectionCoeff, double& diffusionCoeff);

    //@}

  protected:

    //! @name AddNewNeighbours - these functions are collected in AddNewNeighbours for faster search
    //@{

    //! set the neighbours - particles
    virtual void AddNewNeighbours_p(
        const ParticleSPH& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p);

    //! set the neighbours - boundary particles
    virtual void AddNewNeighbours_bp(
        const ParticleSPH& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p);

    //! set the neighbours - free-surface particles
    virtual void AddNewNeighbours_fsp(
        const ParticleSPH& particle_i,
        const std::list<DRT::Node*>& neighboursLinf_p,
        const bool mark_fs_indirect=false);

    //@}

    //! @name interaction methods
    //@{

    //! construct the specific coefficient
    virtual double Inter_generalCoeff_ij(
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate pressure
    virtual void Inter_pressure(
        LINALG::Matrix<3,1>* accn_ij,
        const double& generalCoeff_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate background pressure
    virtual void Inter_backgroundPressure(
        LINALG::Matrix<3,1>* trvl_accn_ij,
        const double& generalCoeff_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate xsph contribution
    virtual void Inter_xsph(
        LINALG::Matrix<3,1>* trvl_accn_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate additional term divA due to transport velocity
    virtual void Inter_transportVelocity_divA(
        LINALG::Matrix<3,1>* accn_ij,
        LINALG::Matrix<3,1>* accn_ij_A,
        const double& generalCoeff_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate laminar viscosity
    virtual void Inter_laminarViscosity(
        LINALG::Matrix<3,1>* accn_ij,
        const double& generalCoeff_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate artificial viscosity
    virtual void Inter_artificialViscosity(
        LINALG::Matrix<3,1>* accn_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
        const InterDataPvP& interData_ij
        );

    //! evaluate surface tension
    virtual void Inter_surfaceTension(
        LINALG::Matrix<3,1>* accn_ij,
        const ParticleSPH& particle_i,
        const ParticleSPH& particle_j,
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

    /// temporary storage of particleSPH data
    std::vector<ParticleSPH> colParticles_;

    /// neighbours - particles
    std::vector<boost::unordered_map<int, InterDataPvP> > neighbours_p_;

    /// weight function
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    /// particle interaction variant
    const bool interactionVariant2_;

    /// wall interaction type
    const INPAR::PARTICLE::WallInteractionType wallInteractionType_;

    /// free surface type
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType_;

    /// temporary storage of wallSPH data
    WallSPHData wallSPHData_;

    /// neighbours - boundary particles
    std::map<int,std::map<int, InterDataPvP> > neighbours_bp_;

    /// neighbours - free surface particles
    std::map<int,std::map<int, InterDataPvP> > neighbours_fsp_;

    /// standard map that contains pointers to all boundary particles of column map
    std::map<int,ParticleSPH* > boundaryparticles_;

    /// standard map that contains pointers to all free-surface particles of column map
    std::map<int,ParticleSPH* > freesurfaceparticles_;

    /// number of space dimensions required for weight function definition
    const INPAR::PARTICLE::WeightFunctionDim WF_DIM_;

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

    /// type of surface tension model
    const INPAR::PARTICLE::SurfaceTensionType surfaceTensionType_;

    /// surface tension between fluid and void phase
    double surfaceTension_;

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
#endif /* SRC_DRT_PARTICLE_PARTICLE_SPH_INTERACTION_H_ */
