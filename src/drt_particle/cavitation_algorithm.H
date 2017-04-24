/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.H

\brief Algorithm to control cavitation simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 11/12 |
 *----------------------------------------------------------------------*/
#ifndef CAVITATION_ALGORITHM_H
#define CAVITATION_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "../drt_inpar/inpar_cavitation.H"
#include "../drt_particle/binning_strategy_utils.H"
#include "../drt_lib/drt_element.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>

#define INLINED_ELE_EVAL

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class Fluid;
}
namespace DRT
{
  class Element;

  namespace MESHFREE
  {
    class MeshfreeMultiBin;
  }
}
namespace LINALG
{
  class SparseMatrix;
}

/*----------------------------------------------------------------------*
 | cavitation algorithm                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
/// CAVITATION: Cavitation problem
namespace CAVITATION
{
  // forward declaration
  class BubbleSource;

  class Algorithm : public PARTICLE::Algorithm
  {
  public:

    Algorithm(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    virtual ~Algorithm(){};

    /// do initialization of problem
    virtual void InitCavitation();

    /// read restart
    virtual void ReadRestart(
      int restart  ///< time step from which to restart
      );

    /// compute norms for normalization during conv check
    virtual void ComputePressAndRadiusNorm(double& pressnorm_L2, double& radiusnorm_L2);

    /// prepare relaxation (save data in the beginning)
    virtual void PrepareRelaxation(const int outeriter);

    /// convergence check and Aitken relaxation
    virtual void ConvergenceCheckAndRelaxation(
      const int outeriter,
      const double pressnorm_L2,
      const double radiusnorm_L2,
      bool& converged,
      bool& particlereset
      );

    /// save bubble data for possible time step repetition
    virtual void SaveParticleData();

    /// reset bubble data when time step is repeated
    virtual void ResetParticleData();

    /// setup of the system
    virtual void SetupSystem();

    /// general time loop for cavitation problem
    virtual void Timeloop();

    /// time loop for cavitation problem with sequential staggering (weak coupling)
    virtual void TimeloopSequStaggered();

    /// time loop for cavitation problem with iterative staggering (strong coupling)
    virtual void TimeloopIterStaggered();

    /// test results (if necessary)
    virtual void TestResults(
      const Epetra_Comm& comm  ///< local epetra communicator
      );


  protected:

    /// initialize bubble pressure
    virtual void InitBubblePressure();

    // initialize bubble velocity from fluid velocity
    virtual void InitBubbleVelFromFluidVel();

    /// compute pressure at bubble position
    virtual bool ComputePressureAtBubblePosition(
      DRT::Node* currparticle,
      const LINALG::Matrix<3,1>& particleposition,
      Epetra_SerialDenseMatrix& elemat1,
      Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1,
      Epetra_SerialDenseVector& elevec2,
      Epetra_SerialDenseVector& elevec3
      );

    /// compute pressure at bubble position
    template <DRT::Element::DiscretizationType distype>
    void ComputePressureAtBubblePositionT(
      const DRT::Element* targetfluidele,
      const LINALG::Matrix<3,1>& elecoord,
      Epetra_SerialDenseMatrix& elemat1,
      Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1,
      Epetra_SerialDenseVector& elevec2,
      Epetra_SerialDenseVector& elevec3
      );

    /// compute velocity at bubble position
    virtual bool ComputeVelocityAtBubblePosition(
      DRT::Node* currparticle,
      LINALG::Matrix<3,1>& particleposition,
      Epetra_SerialDenseMatrix& elemat1,
      Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1,
      Epetra_SerialDenseVector& elevec2,
      Epetra_SerialDenseVector& elevec3
      );

    /// prepare time step
    virtual void PrepareTimeStep();

    /// set coupling states for force computation
    virtual void SetCouplingStates(
      Teuchos::ParameterList& p,
      bool init = false
      );

    /// solve the current particle time step
    virtual void Integrate(bool& particlereset);

    /// apply fluid fraction to fluid field
    virtual void SetFluidFraction();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles(bool init = false);

    /// get data from underlying fluid element
    DRT::Element* GetUnderlyingFluidEleData(
      DRT::Node* currparticle,
      LINALG::Matrix<3,1>& particleposition,
      LINALG::Matrix<3,1>& elecoord,
      const Teuchos::ParameterList& p,
      Epetra_SerialDenseVector& elevector1,
      Epetra_SerialDenseVector& elevector2,
      Epetra_SerialDenseVector& elevector3,
      Epetra_SerialDenseVector& elevector4,
      Epetra_SerialDenseVector& elevector5
      );

    /// get data from underlying fluid element
    template <DRT::Element::DiscretizationType distype>
    void GetUnderlyingFluidEleDataT(
      const LINALG::Matrix<3,1>& particleposition,
      const DRT::Element* targetfluidele,
      LINALG::Matrix<3,1>& elecoord,
      const Teuchos::ParameterList& p,
      Epetra_SerialDenseVector& elevector1,
      Epetra_SerialDenseVector& elevector2,
      Epetra_SerialDenseVector& elevector3,
      Epetra_SerialDenseVector& elevector4,
      Epetra_SerialDenseVector& elevector5
      );

    /// calculate radius based on RP equation
    //! For details see Semesterarbeit Christoph Schmidt
    virtual void ComputeRadius();

    /// time integration of bubble radius
    virtual void IntegrateRadius(
      const int bubbleid,
      double& subtimek,
      const double pvapor,
      const double pg0i,
      double& dtsub,
      const double rho_l,
      const double mu_l,
      const double gamma,
      const double fluid_timenp,
      const double fluid_timen,
      const double pfluidnp,
      const double pfluidn,
      const double r_bub_0,
      const double r_bub_k,
      double& r_bub_kp,
      double& r_dot_bub_k,
      const bool allowadaption,
      const double dt_min,
      const double dt_max
      );

    /// solve Rayleigh-Plesset equation
    virtual double f_R_dot(
      const double r_bub,
      const double r_dot_bub,
      const double rho_l,
      const double sigma,
      const double mu_l,
      const double pbk,
      const double pambientk
      );

    /// insert particles into domain
    virtual void ParticleInflow();

    /// update the current time step
    virtual void Update(const bool converged);

    /// output particle time step
    virtual void Output(bool forced_writerestart = false);

    /// update state vectors which do not live in the time integrator to new layout
    virtual void UpdateStates();

    /// setup of bubble sources in each bin
    virtual void BuildBubbleInflowCondition();

    /// setup ghosting of bins, particles & underlying fluid
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap,  ///< rowmap of bins
      std::map<int, std::set<int> >& rowfluideles,   ///< map from bin id to row fluid elements in this bin
      Teuchos::RCP<Epetra_Map> fluidelecolmapold   ///< element column map before extended ghosting
      );

    /// some output to screen about subcycling
    virtual void SubcyclingInfoToScreen();

    /// build connectivity from fluid elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// get adjacent bins to corner, where ijk is in 1st octant
    virtual std::vector<int> AdjacentBinstoCorner(
      int* ijk   ///< ijk is 1st octant
      );

    //! @name methods for computing the void fraction of the dispersed phase in the fluid
    //! For details see Bachelor thesis Marcel Menner
    //@{

    /// fluid fraction calculation of dispersed phase in fluid
    virtual void CalculateFluidFraction(
      Teuchos::RCP<const Epetra_Vector> particleradius
      );

    /// use Gaussian integration for solving volume integral to obtain void volume
    virtual void DoGaussianIntegrationFluidFrac(
      DRT::Element*              ele,
      LINALG::Matrix<3,1>&       particleposition,
      const double               influence,
      double&                    vol_ele,
      const bool                 pbcdetected,
      const LINALG::Matrix<3,1>& pbceleoffset
      );

    /// use analytical integration for solving volume integral to obtain void volume
    virtual void DoAnalyticalIntegrationFluidFrac(
      DRT::Element*               ele,
      const int                   bubbleid,
      const LINALG::Matrix<3,1>&  particleposition,
      const double                influence,
      double&                     vol_ele,
      bool&                       surfaceoverlap,
      const bool                  bcpbcdetected,
      const LINALG::Matrix<3,1>&  pbceleoffset
      );

    /// compute integration points for analytic integration
    virtual int EvaluateSurface(
      std::vector<LINALG::Matrix<3,1> >& surfacenodes,
      const LINALG::Matrix<3,1>&       n,
      const LINALG::Matrix<3,1>&       centerele ,
      const LINALG::Matrix<3,1>&       particleposition,
      const double                     influence,
      double&                          vol_ele,
      bool&                            surfaceoverlap,
      const int                        eleid,
      const int                        isurface
      );

    /// find penetration points of bubble surfaces (+x and -x)
    virtual void GetPenetrationPointsOfXSurfaces(
      const LINALG::Matrix<3,1>&         n,
      std::vector<LINALG::Matrix<3,1> >& surfacenodes,
      const double                       influence,
      const LINALG::Matrix<3,1>&         particleposition,
      std::vector<LINALG::Matrix<3,1> >& bubbleX,
      std::vector<LINALG::Matrix<3,1> >& bubble_X,
      bool&                              surfaceoverlap,
      const int                          eleid
      );

    /// calculation of bubble corner points
    virtual void CalculateBubbleCornerPoint(
      const LINALG::Matrix<3,1>&    n,
      const LINALG::Matrix<3,1>&    centersurface,
      int                           bubblesurface1,
      int                           bubblesurface2,
      double                        *bubblesurface,
      LINALG::Matrix<3,1>&          bubblecorner
      );

    /// check whether a given point is inside of a surface
    virtual bool CheckPointInSurface(
      std::vector<LINALG::Matrix<3,1> >& surfacenodes,
      const LINALG::Matrix<3,1>&         centersurface,
      const LINALG::Matrix<3,1>&         centerele,
      const LINALG::Matrix<3,1>&         pointtocheck
      );

    /// build convex hull of points in y-z plane
    virtual void BuildConvexHull(
      std::vector<LINALG::Matrix<3,1> >& surfacenodes
      );

    /// evaluate line integral for analytic integration (constant polynomial)
    virtual void EvaluateTwoPointsConstPoly(
      const LINALG::Matrix<3,1>&    n,
      const LINALG::Matrix<3,1>&    CenterRing,
      const LINALG::Matrix<3,1>&    point1,
      const LINALG::Matrix<3,1>&    point2,
      double&                       vol_ele
      );

    /// evaluate line integral for analytic integration (quadratic polynomial)
    virtual void EvaluateTwoPointsQuadraticPoly(
      const LINALG::Matrix<3,1>&    n,
      const LINALG::Matrix<3,1>&    centerringintgral,
      LINALG::Matrix<3,1>&          point1,
      LINALG::Matrix<3,1>&          point2,
      const LINALG::Matrix<3,1>&    particleposition,
      double&                       vol_ele,
      const double                  influence
      );

    /// evaluate line integral for analytic integration (quartic polynomial)
    virtual void EvaluateTwoPointsQuarticPoly(
      const LINALG::Matrix<3,1>& n,
      const LINALG::Matrix<3,1>& centerringintgral,
      LINALG::Matrix<3,1>& point1,
      LINALG::Matrix<3,1>& point2,
      const LINALG::Matrix<3,1>& particleposition,
      double& vol_ele,
      const double influence
      );

    //@}

    /// test for overlap of the XAABB of bubble and element and compute pbc offset if necessary
    virtual bool XAABBoverlap(
      DRT::Element* ele,
      const double influence,
      const LINALG::Matrix<3,1>& particleposition,
      const bool havepbc,
      bool& pbcdetected,
      LINALG::Matrix<3,1>& pbceleoffset
      );

    /// assign bubble volume to the underlying fluid element
    virtual void AssignSmallBubbles(
      const double bubblevol,
      const LINALG::Matrix<3,1>& particleposition,
      const std::vector<int>& insideeles,
      Teuchos::RCP<Epetra_FEVector> void_volumes,
      const std::vector<bool>& pbcdetected,
      const std::vector<LINALG::Matrix<3,1> >& pbceleoffset
      );

    /// get underlying element as well as position in element space
    virtual DRT::Element* GetEleCoordinatesFromPosition(
      const DRT::Node* currparticle,
      const LINALG::Matrix<3,1>& myposition,
      DRT::MESHFREE::MeshfreeMultiBin* currbin,
      LINALG::Matrix<3,1>& elecoord,
      const bool approxelecoordsinit
      );

    /// unwarp surfaces of fluid element for analytical integration
    virtual void UnwarpElement(
      const DRT::Element* ele,
      const std::vector<Teuchos::RCP<DRT::Element> >& surfaces,
      const int numsurfacenodes,
      std::map<int, LINALG::Matrix<3,1> >& currentpositions,
      std::vector<LINALG::Matrix<3,1> >& normals
      );

    /// compute superconvergent patch recovery for fluid fraction
    virtual void ComputePatchRecoveredFluidFraction(
      Teuchos::RCP<const Epetra_MultiVector> fluidfraction
      );

    /// L2 projection of a single value per element to nodal value distribution
    virtual void ComputeL2ProjectedFluidFraction(
      Teuchos::RCP<const Epetra_MultiVector>
      );

    /// solve for a node based linear system
    virtual Teuchos::RCP<Epetra_Vector> SolveOnNodeBasedVector(
      const int solvernumber,
      Teuchos::RCP<Epetra_Map> noderowmap,
      Teuchos::RCP<LINALG::SparseMatrix> matrix,
      Teuchos::RCP<Epetra_Vector> rhs
      );

    /// print cut situation for analytical void fraction computation to GMSH
    virtual void PrintBubbleAndFluidEleToGMSH(
      DRT::Element* ele,
      const int bubbleid,
      const LINALG::Matrix<3,1>& particleposition,
      const double influence
      );

    /// print integration lines for analytical void fraction computation to GMSH
    virtual void PrintIntegrationLinesToGMSH(
      const std::vector<LINALG::Matrix<3,1> >& integrationpoints,
      const int isurface
      );

  private:
    const int dim_;   ///< spatial dimension of problem

    //! @name cavitation related stuff
    //@{

    const INPAR::CAVITATION::CouplingStrategyOverFields coupalgo_;   ///< coupling strategy between fluid and particle

    const INPAR::CAVITATION::VoidFractionCalculation void_frac_strategy_;   ///< strategy for void fraction calculation

    const INPAR::CAVITATION::FluidFracReconstructionMethod fluidfrac_reconstr_;   ///< type of fluid fraction reconstruction

    const int gauss_rule_per_dir_;   ///< number of gauss points in each direction for void fraction computation

    const bool approxelecoordsinit_;   ///< bool whether a fast and approximate initial guess is used for computing natural coords

    const bool simplebubbleforce_;   ///< bool whether a simplified force computation is applied

    const int timestepsizeratio_;   ///< ratio between fluid and particle time step in case of subcycling

    int restartparticles_;   ///< restart step for particle subproblem

    bool inflowradiusblending_;   ///< bool whether radius of inflow particles is blended over several time steps

    int blendingsteps_;   ///< number of steps that is used to blend bubble radius from initially small to final size

    bool initbubblevelfromfluid_;   ///< bool whether initial bubble velocity is computed from the underlying fluid element

    const double influencescaling_;   ///< scaling factor for influence of bubble in void frac computation

    //@}

    //! @name fluid related stuff
    //@{

    BINSTRATEGY::UTILS::BinContentType bin_fluidcontent_;

    Teuchos::RCP<DRT::Discretization> fluiddis_;     ///< underlying fluid discretization

    Teuchos::RCP<ADAPTER::Fluid> fluid_;     ///< underlying fluid time integrator

    std::map<int, std::list<Teuchos::RCP<BubbleSource> > > bubble_source_;   ///< store particle sources for each bin

    std::set<int> latestinflowbubbles_;   ///< store latest inflow particle ids

    std::set<int> inflowfluiddofs_;  /// variable to store all dofs of fluid nodes close to the inflow

    Teuchos::RCP<Epetra_Vector> ele_volume_;   ///< element row map containing the volume of each element

    Teuchos::RCP<Epetra_Vector> fluidfracn_;   ///< fluid fraction in pressure dof of dofrowmap at t_n

    Teuchos::RCP<Epetra_Vector> fluidfracnp_;   ///< fluid fraction in pressure dof of dofrowmap at t_n+1

    bool fluidfrac_relevant;   ///< bool to check whether fluid fraction is considered in the problem type

    Teuchos::RCP<Epetra_Vector> velcoln_;   ///< velocity in col layout at t_n

    Teuchos::RCP<Epetra_Vector> velcolnp_;   ///< velocity in col layout at t_n+1

    Teuchos::RCP<Epetra_Vector> acccolnp_;   ///< acceleration in col layout at t_n+1

    Teuchos::RCP<Epetra_Vector> velcol_interpol_;   ///< velocity in col layout at t_{n+theta}

    Teuchos::RCP<Epetra_MultiVector> velgradcoln_;   ///< velocity gradient in col layout at t_n

    Teuchos::RCP<Epetra_MultiVector> velgradcolnp_;   ///< velocity gradient in col layout at t_n+1

    Teuchos::RCP<Epetra_MultiVector> velgradcol_interpol_;   ///< velocity gradient in col layout at t_{n+theta}

    //@}

    //! @name Rayleigh-Plesset specific stuff for bubble radius variation
    //@{

    const bool computeradiusRPbased_;   ///< bool whether bubble radius is adapted based on the Rayleigh-Plesset equation

    Teuchos::RCP<Epetra_Vector> dtsub_;   ///< individual time step size for each bubble during radius adaption

    Teuchos::RCP<Epetra_Vector> pg0_;   ///< partial pressure of noncondensable gas at bubble initialization time

    //@}

    //! @name performance related stuff
    //@{

    int count_;   ///< counter for triggering particle transfer

    struct UnderlyingEle{
      DRT::Element* ele;
      LINALG::Matrix<3,1> elecoord;
    };

    std::vector<UnderlyingEle> underlyingelecache_;   ///< cache for underlying fluid element

    std::vector<LINALG::SerialDenseMatrix> xyze_cache_;   ///< prefetched coordinates of fluid element
    std::vector<std::vector<int> > lm_cache_;   ///< prefetched lm vectors of fluid elements
    std::map<int, Epetra_SerialDenseVector> evelgrad_cache_;   ///< cached evel grad of fluid elements
    std::map<int, std::vector<std::vector<double> > > evelacc_cache_;   ///< cached evel and eacc of fluid elements

    //@}

    //! @name Aitken relaxation

    Teuchos::RCP<Epetra_Vector>  del_;   ///< difference of last solution and relaxed solution which was used to compute last solution
                                         ///< del = r^{i+1}_{n+1} = d^{i+1}_{n+1} - d^{i,relaxed}_{n+1}

    Teuchos::RCP<Epetra_Vector> delhist_;   ///< difference of difference of last two pair of solutions
                                            ///< delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )

    double mu_;   ///< Aitken factor

    Teuchos::RCP<Epetra_Vector> radius_i_;   ///< radius of last outer iteration for dynamic relaxation

    Teuchos::RCP<Epetra_Vector> couplingradius_;   ///< radius at current outer iteration for dynamic relaxation

    Teuchos::RCP<Epetra_Vector> pressnp_;   ///< fluid pressure for convergence check

    const int itmax_;   ///<  maximum iteration steps

    const double ittol_;   ///< convergence tolerance

    /// safe data for repeating time steps
    int storecount_;
    double storetime_;
    int storestep_;
    Teuchos::RCP<const Epetra_Vector> storedis_;
    Teuchos::RCP<const Epetra_Vector> storevel_;
    Teuchos::RCP<const Epetra_Vector> storeacc_;
    Teuchos::RCP<const Epetra_Vector> storeang_vel_;
    Teuchos::RCP<const Epetra_Vector> storeang_acc_;
    Teuchos::RCP<const Epetra_Vector> storerad_;
    Teuchos::RCP<const Epetra_Vector> storeraddot_;
    Teuchos::RCP<const Epetra_Vector> storemass_;
    Teuchos::RCP<const Epetra_Vector> storeinertia_;
    Teuchos::RCP<const Epetra_Vector> storedtsub_;
    Teuchos::RCP<const Epetra_Vector> storerad0_;
    Teuchos::RCP<const Epetra_Vector> storepg0_;
    std::list<int> deletedparticlesduringsubcycling_;

    //@}

    /// gmsh output
    std::ofstream gmshfilecontent_;


  }; // Algorithm


  class BubbleSource
  {
  public:

    BubbleSource(
      const int bubbleinflowid,   ///< unique id over all inflow conditions
      std::vector<double> inflow_position,  ///< position of particle source
      std::vector<double> inflow_vel,   ///< inflow velocity
      const int inflow_vel_curve,   ///< inflow velocity can be superposed with a time curve
      const double inflow_radius,   ///< initial radius of particle
      const double inflow_freq,   ///< frequency of entering particles
      const double timedelay,   ///< time delay for inflowing particles
      const double stopinflowtime   ///< time when inflow should stop
      );

    // internal variables
    const int inflowid_;
    std::vector<double> inflow_position_;
    std::vector<double> inflow_vel_;
    const int inflow_vel_curve_;
    const double inflow_radius_;
    const double inflow_freq_;
    const double timedelay_;
    const double stopinflowtime_;

  }; // BubbleSource

} // namespace CAVITATION


/*----------------------------------------------------------------------*/
#endif  // CAVITATION_ALGORITHM_H
