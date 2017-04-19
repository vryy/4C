/*!-----------------------------------------------------------------------------------------------*
\file levelset_algorithm.H

\brief level-set algorithm

\level 2

\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

 *------------------------------------------------------------------------------------------------*/

#ifndef LEVELSET_ALGORITHM_H
#define LEVELSET_ALGORITHM_H

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_inpar/inpar_levelset.H"
#include "../drt_geometry/geo_utils.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_fixedsizematrix.H"

#define USE_PHIN_FOR_VEL // TODO


namespace PARTICLE
{
  class ScatraParticleCoupling;
}

namespace LINALG
{
  class SerialDenseMatrix;
}


namespace SCATRA
{

  class LevelSetAlgorithm : public virtual ScaTraTimIntImpl
  {
  public:

    /// Standard Constructor
    LevelSetAlgorithm(
        Teuchos::RCP<DRT::Discretization>        dis,
        Teuchos::RCP<LINALG::Solver>             solver,
        Teuchos::RCP<Teuchos::ParameterList>     params,
        Teuchos::RCP<Teuchos::ParameterList>     sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList>     extraparams,
        Teuchos::RCP<IO::DiscretizationWriter>   output);

    /// virtual destructor
    virtual ~LevelSetAlgorithm();

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// initialize level-set algorithm
    virtual void Init();

    /// setup level-set algorithm
    virtual void Setup();

    /// time loop
    void TimeLoop();

    /// output solution and restart data to file
    void Output(const int num =0);

    /// read restart data
    virtual void ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null) = 0;

    //! set the velocity field (zero or field by function) (pure level-set problems)
    void SetVelocityField(bool init=false);

    /// set convective velocity field (+ pressure and acceleration field as
    /// well as fine-scale velocity field, if required) (function for coupled fluid problems)
    void SetVelocityField(
      Teuchos::RCP<const Epetra_Vector> convvel,
      Teuchos::RCP<const Epetra_Vector> acc,
      Teuchos::RCP<const Epetra_Vector> vel,
      Teuchos::RCP<const Epetra_Vector> fsvel,
      bool setpressure=false,
      bool init=false);

    /// access routine for nodal curvature
    Teuchos::RCP<Epetra_Vector> GetNodalCurvature(const Teuchos::RCP<const Epetra_Vector> phi,
                                                  const Teuchos::RCP<const Epetra_MultiVector> gradphi);


    // output position of center of mass assuming a smoothed interfaces
    void MassCenterUsingSmoothing();

    /// redistribute the scatra discretization and vectors according to nodegraph
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> & nodegraph);

    /// return velocity at intermediate time n+theta
    const Teuchos::RCP<Epetra_Vector> ConVelTheta(double theta);

    /// perform result test
    void TestResults();

    /// set time and step value
    void SetTimeStep(const double time, const int step);

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// solve level-set equation
    void Solve();

    /// calculate error compared to analytical solution
    void EvaluateErrorComparedToAnalyticalSol();

  protected:

    virtual void GetInitialVolumeOfMinusDomain(
        const Teuchos::RCP<const Epetra_Vector>& phinp,
        const Teuchos::RCP<const DRT::Discretization> & scatradis,
        double& volumedomainminus) const;

    /// update state vectors
    virtual void UpdateState() = 0;

    // -----------------------------------------------------------------
    // hybrid method including particles
    // -----------------------------------------------------------------

    /// correct zero level-set by particles
    void ParticleCorrection();

    // -----------------------------------------------------------------
    // reinitialization
    // -----------------------------------------------------------------

    /// reinitialize level-set
    virtual void Reinitialization();

    /** \brief calculation of nodal velocity field via L2-projection for reinitialization
     *
     * (helper function for ReinitEq()) */
    virtual void CalcNodeBasedReinitVel();

    /// execute the elliptic reinitialization procedure
    void ReinitializeWithEllipticEquation();

    /// set element parameters for reinitialization equation
    void SetReinitializationElementParameters(bool calcinitialtimederivative = false) const;

    /** \brief access nodal gradient-based values for reinitialization
     *
     * (ReinitEq() only; Sussman and Elliptic) */
    inline Teuchos::RCP<Epetra_MultiVector> & NodalGradBasedValue()
    { return nb_grad_val_; }
    inline Teuchos::RCP<const Epetra_MultiVector> NodalGradBasedValue() const
    { return nb_grad_val_; }

    /// access the interval for reinitialization (every 'reinitinterval_' time steps)
    inline const int & ReinitInterval() const
    { return reinitinterval_; }

    /// use projection for grad phi and related quantities
    inline const bool & UseProjection() const
    { return projection_; }

    // -----------------------------------------------------------------
    // Reconstructing nodal curvature
    // -----------------------------------------------------------------
    void ReconstructedNodalCurvature(
      Teuchos::RCP<Epetra_Vector> curvature,
      const Teuchos::RCP<const Epetra_Vector> phi,
      const Teuchos::RCP<const Epetra_MultiVector> gradphi);

    // -----------------------------------------------------------------
    // members
    // -----------------------------------------------------------------

    /// the parameter list for level-set problems
    Teuchos::RCP<Teuchos::ParameterList> levelsetparams_;

    /// pointer to particle algorithm
    Teuchos::RCP<PARTICLE::ScatraParticleCoupling> particle_;

    /// options for reinitialization of G-function;
    INPAR::SCATRA::ReInitialAction reinitaction_;

    /// convective velocity at time n
    Teuchos::RCP<const Epetra_Vector> conveln_;

    /// flag to switch between standard integration and sub-time loop for reinitialization
    bool switchreinit_;

    /// maximal number of pseudo time steps (ReinitEq() only)
    int pseudostepmax_;

    /// pseudo time step counter (ReinitEq() only)
    int pseudostep_;

    /// pseudo time step length (ReinitEq() only)
    double dtau_;

    /// pseudo theata (ReinitEq() only)
    double thetareinit_;

  private:

    /// add parameters depending of the problem, i.e., loma, level-set, ...
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params);

    /// manipulate velocity field away from the interface
    void ManipulateFluidFieldForGfunc();

    /// modification of convective velocity at contact points
    void ApplyContactPointBoundaryCondition();

    // -----------------------------------------------------------------
    // reinitialization
    // -----------------------------------------------------------------

    /// algebraic reinitialization via solution of equation to steady state
    void ReinitEq();

    /// set time parameters for reinitialization equation
    void SetReinitializationElementTimeParameters();

    /// preparations to solve reinitialization equation within existing framework (helper function for ReinitEq())
    void PrepareTimeLoopReinit();

    /// time loop for reinitialization equation (helper function for ReinitEq())
    void TimeLoopReinit();

    /// clean necessary modifications to solve reinitialization equation within existing framework (helper function for ReinitEq())
    void FinishTimeLoopReinit();

    /// setup the variables to do a new reinitialization time step (helper function for ReinitEq())
    void PrepareTimeStepReinit();

    /// nonlinear solver for reinitialization equation (helper function for ReinitEq())
    void SolveReinit();

    /// correction step according to Sussman & Fatemi 1999 (helper function for ReinitEq())
    void CorrectionReinit();

    /// convergence check for reinit equation according to Sussman et al. 1994 (helper function for ReinitEq())
    bool ConvergenceCheckReinit();

    /// update phi within the reinitialization loop
    virtual void UpdateReinit() = 0;

    /// geometric reinitialization via computation of distance of node to interface
    void ReinitGeo(const std::map<int,GEO::BoundaryIntCells >& interface);

    /// compute normal vector of interface patch (helper function for ReinitGeo())
    void ComputeNormalVectorToInterface(
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        LINALG::Matrix<3,1>&             normal);

    /// compute distance to vertex of patch (helper function for ReinitGeo())
    void ComputeDistanceToPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        double&                          vertexdist);

    /// compute distance to edge of patch (helper function for ReinitGeo())
    void ComputeDistanceToEdge(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        double&                          edgedist);

    /// find a facing interface patch by projection of node into boundary cell space (helper function for ReinitGeo())
    void FindFacingPatchProjCellSpace(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        bool&                            facenode,
        double&                          patchdist);

    /// compares the second entry of a pair<int,double>. To be passed to the sorting algo (helper function for ReinitGeo())
    static bool MyComparePairs (const std::pair<int,double>& first, const std::pair<int,double>& second)
    {
      if (fabs(first.second) < fabs(second.second))
        return true;
      else
        return false;
    };

    /// project node into the boundary cell space (2D) (helper function for ReinitGeo())
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ProjectNodeOnPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        LINALG::Matrix<2,1>&             eta,
        double&                          alpha);

    /// correct the volume of the minus domain after reinitialization
    void CorrectVolume();

    /// elliptic reinitialization
    void ReinitElliptic(
        std::map<int,GEO::BoundaryIntCells >& interface);

    // -----------------------------------------------------------------
    // additional post-processing and evaluation methods
    // -----------------------------------------------------------------

    virtual void OutputOfLevelSetSpecificValues();

    // check for mass conservation before and after reinitialization as well as
    // at the end of the time step
    void MassConservationCheck(
      const double actvolminus,
      const bool writetofile = false);

    // reconstruction of interface and output of domains
    void CaptureInterface(
      std::map<int,GEO::BoundaryIntCells >& interface,
      const bool writetofile = false);

    // -----------------------------------------------------------------
    // members
    // -----------------------------------------------------------------

    // --------------
    // members related to reinitialization and mass conservation

    /// initial volume of minus domain
    double initvolminus_;

    /// phinp before reinitialization (ReinitEq() only)
    Teuchos::RCP<Epetra_Vector> initialphireinit_;

    /// nodal gradient-based values for reinitialization (ReinitEq() only; Sussman and Elliptic)
    Teuchos::RCP<Epetra_MultiVector> nb_grad_val_;

    /// interval for reinitialization (every 'reinitinterval_' time steps)
    int reinitinterval_;

    /// switch for reinitialization only within a band around the interface (ReinitGeo() only)
    bool reinitband_;

    /// band width for reinitialization (maximum level-set value)
    double reinitbandwidth_;

    /// flag to activate corrector step (ReinitEq() only)
    bool reinitcorrector_;

    /// evaluation of velocity field for reinitialization (ReinitEq() only)
    INPAR::SCATRA::VelReinit useprojectedreinitvel_;

    /// 2D flag
    INPAR::SCATRA::LSDim lsdim_;

    /// use projection for grad phi and related quantities
    bool projection_;

    //TODO:
//    /// vector containing denominator of penalty parameter for each element (ReinitEq() only)
//    Teuchos::RCP<Epetra_Vector> lambda_ele_denominator_;
//
//    /// vector containing smoothed haevyside function for each dof (ReinitEq() only)
//    Teuchos::RCP<Epetra_Vector> node_deriv_smoothfunct_;

    /// tolerance for convergence check according to Sussman et al. 1994 (turned off negative) (ReinitEq() only)
    double reinit_tol_;

    /// flag to correct volume after reinitialization
    bool reinitvolcorrection_;

    /// interface for elliptic reinitialization
    Teuchos::RCP<std::map<int,GEO::BoundaryIntCells > > interface_eleq_;

    // --------------
    // members related to transport and velocity fields from Navier-Stokes

    /// flag for replacing the velocity of nodes at a certain distance of the interface by an approximate interface velocity
    bool extract_interface_vel_;

    /// number of element layers around the interface with unmodified velocity
    int convel_layers_;

    /// flag for modification of convective velocity at contact points
    bool cpbc_;

  };

}

#endif // #ifndef LEVELSET_ALGORITHM_H
