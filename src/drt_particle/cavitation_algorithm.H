/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.H

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

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

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>

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

    /// setup of the system
    virtual void SetupSystem();

    /// outer time loop for cavitation problem
    virtual void Timeloop();

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
      LINALG::Matrix<3,1>& particleposition,
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

    /// solve the current particle time step
    virtual void Integrate();

    /// apply fluid fraction to fluid field
    virtual void SetFluidFraction();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles();

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
    virtual void Update();

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
    virtual void CalculateFluidFraction();

    /// use Gaussian integration for solving volume integral to obtain void volume
    virtual void DoGaussianIntegrationFluidFrac(
      DRT::Element*         ele,
      LINALG::Matrix<3,1>&  particleposition,
      const double          influence,
      double&               vol_ele
      );

    /// use analytical integration for solving volume integral to obtain void volume
    virtual void DoAnalyticalIntegrationFluidFrac(
      DRT::Element*         ele,
      LINALG::Matrix<3,1>&  particleposition,
      const double          influence,
      double&               vol_ele,
      bool&                 surfaceoverlap
      );

    /// compute integration points for analytic integration
    virtual void EvaluateSurface(
      std::vector<LINALG::Matrix<3,1> >& surfacenodes,
      const LINALG::Matrix<3,1>&       n,
      const LINALG::Matrix<3,1>&       centerele ,
      LINALG::Matrix<3,1>&             particleposition,
      const double                     influence,
      double&                          vol_ele,
      bool&                            surfaceoverlap
      );

    /// find penetration points of bubble surfaces (+x and -x)
    virtual void GetPenetrationPointsOfXSurfaces(
      const LINALG::Matrix<3,1>&         n,
      std::vector<LINALG::Matrix<3,1> >& surfacenodes,
      const double                       influence,
      const LINALG::Matrix<3,1>&         particleposition,
      std::vector<LINALG::Matrix<3,1> >& bubbleX,
      std::vector<LINALG::Matrix<3,1> >& bubble_X,
      bool&                              surfaceoverlap
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

    //@}

    /// test for overlap of the XAABB of bubble and element
    virtual bool XAABBoverlap(
      DRT::Element* ele,
      const double influence,
      LINALG::Matrix<3,1> particleposition
      );

    /// assign bubble volume to the underlying fluid element
    virtual void AssignSmallBubbles(
      const double bubblevol,
      const LINALG::Matrix<3,1> particleposition,
      const std::vector<int> insideeles,
      Teuchos::RCP<Epetra_FEVector> void_volumes
      );

    /// get underlying element as well as position in element space
    virtual DRT::Element* GetEleCoordinatesFromPosition(
      LINALG::Matrix<3,1>& myposition,
      DRT::MESHFREE::MeshfreeMultiBin* currbin,
      LINALG::Matrix<3,1>& elecoord,
      const bool approxelecoordsinit
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

    //@}

    //! @name fluid related stuff
    //@{

    Teuchos::RCP<DRT::Discretization> fluiddis_;     ///< underlying fluid discretization

    Teuchos::RCP<ADAPTER::Fluid> fluid_;     ///< underlying fluid time integrator

    std::map<int, std::list<Teuchos::RCP<BubbleSource> > > bubble_source_;   ///< store particle sources for each bin

    std::set<int> latestinflowbubbles_;   ///< store latest inflow particle ids

    std::set<int> inflowfluiddofs_;  /// variable to store all dofs of fluid nodes close to the inflow

    Teuchos::RCP<Epetra_Vector> ele_volume_;   ///< element row map containing the volume of each element

    Teuchos::RCP<Epetra_Vector> fluidfracn_;   ///< fluid fraction in pressure dof of dofrowmap at t_n

    Teuchos::RCP<Epetra_Vector> fluidfracnp_;   ///< fluid fraction in pressure dof of dofrowmap at t_n+1

    //@}

    //! @name Rayleigh-Plesset specific stuff for bubble radius variation
    //@{

    const bool computeradiusRPbased_;   ///< bool whether bubble radius is adapted based on the Rayleigh-Plesset equation

    Teuchos::RCP<Epetra_Vector> dtsub_;   ///< individual time step size for each bubble during radius adaption

    Teuchos::RCP<Epetra_Vector> pg0_;   ///< partial pressure of noncondensable gas at bubble initialization time

    //@}

    int count_;   ///< counter for triggering particle transfer

  }; // Algorithm


  class BubbleSource
  {
  public:

    BubbleSource(
      int bubbleinflowid,   ///< unique id over all inflow conditions
      std::vector<double> inflow_position,  ///< position of particle source
      std::vector<double> inflow_vel,   ///< inflow velocity
      int inflow_vel_curve,   ///< inflow velocity can be superposed with a time curve
      double inflow_radius,   ///< initial radius of particle
      double inflow_freq,   ///< frequency of entering particles
      double timedelay
      );

    // internal variables
    int inflowid_;
    std::vector<double> inflow_position_;
    std::vector<double> inflow_vel_;
    int inflow_vel_curve_;
    double inflow_radius_;
    double inflow_freq_;
    double timedelay_;

  }; // BubbleSource

} // namespace CAVITATION


/*----------------------------------------------------------------------*/
#endif  // CAVITATION_ALGORITHM_H
