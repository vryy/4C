/*!----------------------------------------------------------------------
\file combust_fluidimplicitintegration.H
\brief class holding implicit time integration schemes for combustion problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#ifndef COMBUST_FLUIDIMPLICITINTEGRATION_H
#define COMBUST_FLUIDIMPLICITINTEGRATION_H


#include "../drt_fluid/fluid_timint.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_io/io_pstream.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/xfem_utils_combust.H"

#include <Teuchos_ParameterList.hpp>

// forward declarations
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Vector;
class PeriodicBoundaryConditions;

namespace COMBUST
{
  class FlameFront;
  class InterfaceHandleCombust;
}

namespace DRT
{
  class Discretization;
  class DiscretizationCombust;
  class IndependentPBCDofSet;
}

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
  class Solver;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class DofKey;
  class DofManager;
  class ElementAnsatz;
  class ENR;
  class STD;
}

namespace FLD
{
  class TransferTurbulentInflowCondition;

/*!
\brief implicit time integration for combustion problems

\author henke (rasthofer@lnm.mw.tum.de)
*/

// forward declarations
class TurbulenceStatisticManager;

class CombustFluidImplicitTimeInt : public TimInt
{
  friend class TurbulenceStatisticManager;
  friend class CombustFluidResultTest;

public:

  /*!
  \brief constructor

  */
  explicit CombustFluidImplicitTimeInt(
      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
    );

  /*!
  \brief destructor

  */
  virtual ~CombustFluidImplicitTimeInt();

  /// initialization of all members
  void Init();

  /// return type of time integration scheme
  INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

  /// Provide title of time integration scheme
  std::string TimIntTitle(){return MapTimIntEnumToString(TimIntScheme());}

  /// Print information about current time step to screen
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
    {
      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta)
      {
        IO::cout << "TIME: "
                 << std::setw(11) << std::setprecision(4) << std::scientific << time_ << "/"
                 << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_ << "  DT = "
                 << std::setw(11) << std::setprecision(4) << std::scientific << dta_ << "  "
                 << TimIntTitle() << " (theta = "
                 << std::setw(3) << theta_ << ") STEP = "
                 << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_ << IO::endl;
      }
      else{
        IO::cout << "TIME: "
                 << std::setw(11) << std::setprecision(4) << std::scientific << time_ << "/"
                 << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_ << "  DT = "
                 << std::setw(11) << std::setprecision(4) << std::scientific << dta_ << "  "
                 << TimIntTitle() << "  STEP = "
                 << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_ << IO::endl;
      }
    }
  }

  /*!
  \brief solve stationary problem
         Even though it does not make much sense, this member function belongs to the implicit time
         integration scheme for convenience. Of course, there is no time integration for stationary
         problems.
  */
  void SolveStationaryProblem();


  /// check whether FRS iteration is finished
  bool FluidRefSolLoopFinished();

  /// clear time integration data for XFEM problems
  void ClearTimeInt();

  /// prepare a time step
  void PrepareTimeStep();

  /// prepare a nonlinear iteration
  void PrepareSolve();

  /// import the flame front from the combustion algorithm
  void StoreFlameFront(
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      bool UpdateDofSet);

  /// Implement ADAPTER::FLUID
  void ImportFlameFront(const Teuchos::RCP<COMBUST::FlameFront> flamefront, bool UpdateDofSet)
  {
    StoreFlameFront(flamefront, UpdateDofSet);
  }

  /*!
  \brief get the velocity vector based on standard dofs

  \return Teuchos::RCP to a copy of Veln with only standard dofs
   */
  Teuchos::RCP<Epetra_Vector> StdVeln();

  /*!
  \brief get the velocity vector based on standard dofs

  \return Teuchos::RCP to a copy of Velnp with only standard dofs
   */
  Teuchos::RCP<Epetra_Vector> StdVelnp();

  /*!
  \brief get the velocity vector based on standard dofs

  \return Teuchos::RCP to a copy of Velaf with only standard dofs
   */
  Teuchos::RCP<Epetra_Vector> StdVelaf();

  /// get history vector for transfer to scalar transport field
  Teuchos::RCP<const Epetra_Vector> Hist();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or fixpoint iteration
  */
  void Solve();

  Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; };

  /*!
  \brief preparations for Krylov space projection
  */
  void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
  void UpdateKrylovSpaceProjection();
  void CheckMatrixNullspace();

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  <pre>

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    //
    // One-step-Theta: (step>1)
    //
    //  accn_  = (velnp_-veln_) / (Theta * dt) - (1/Theta -1) * accn_
    //  "(n+1)"
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2:           (step>1)
    //
    //               2*dt(n)+dt(n-1)              dt(n)+dt(n-1)
    //  accn_   = --------------------- velnp_ - --------------- veln_
    //             dt(n)*[dt(n)+dt(n-1)]          dt(n)*dt(n-1)
    //
    //                     dt(n)
    //           + ----------------------- velnm_
    //             dt(n-1)*[dt(n)+dt(n-1)]
    //
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2 and  One-step-Theta: (step==1)
    //
    // The given formulas are only valid from the second timestep. In the
    // first step, the acceleration is calculated simply by
    //
    //  accn_  = (velnp_-veln_) / (dt)
    //
    // -------------------------------------------------------------------

  </pre>

  */
  void TimeUpdate();

  /// Implement ADAPTER::Fluid
  void Update() { TimeUpdate(); }

  /*
  \brief provide access to turbulence statistics manager
  */
  Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager() { return turbstatisticsmanager_; }
  /// provide access to the box filter for dynamic Smagorinsky model
  Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter(){return Teuchos::null;}

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  void SetInitialFlowField(const INPAR::FLUID::InitialField initfield, const int startfuncno)
  {
    // This function is called from the Fluid Base Algorithm, but has no effect.
    // The fluid flow field is initialized from the COMBUST::Algorithm.
    return;
  }

  /*!
  \brief set initial flow field, also for analytical test problems
  */
  void SetInitialFlowField(const INPAR::COMBUST::InitialField initfield, const int initfuncno);

  /*!
  \brief set initial enrichment field for analytical test problems
  */
  void SetEnrichmentField(
      const Teuchos::RCP<XFEM::DofManager> dofmanager,
      const Epetra_Map dofrowmap);

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems
  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief calculate error between an analytical solution and the numerical solution of a test problem
  */
  void EvaluateErrorComparedToAnalyticalSol_Nitsche(INPAR::COMBUST::NitscheError& );

  /*!
  \brief calculate stresses
  */
//  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief out of use
  */
  void Integrate();

  /*!
  \brief out of use
  */
  void TimeLoop();

  /*!
  \brief not used in combustion or XFEM -> remove
         build linear system matrix and rhs
  \param vel new guess at velocity and pressure
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);

  /*!
  \brief output solution fields to Gmsh for debugging

  */
  void OutputToGmsh(
      const char* presName,
      const char* velName,
      const int step,
      const double  time
  ) const;

  /*!
  \brief compute flame area and output it to file

  */
  void OutputFlameArea(
      const int     step,
      const double  time
  ) const;

  /*!
  \brief compute amplitude of instability and output it to file

  */
  void OutputInstabAmplitude(
      const int     step,
      const double  time
  ) const;

  /*!
   * \brief plot a vector field to a Gmsh postprocessing file
   *
   * @param vectorfield      vector with solution values
   * @param filestr          basic name of the output file
   * @param name_in_gmsh     text that appears in the Gmsh window over the scale
   * @param plot_to_gnuplot  plot some values to Gnuplot (debugging - will go away sooner or later)
   */
  void PlotVectorFieldToGmsh(
      const Teuchos::RCP<const Epetra_Vector>   vectorfield,
      const std::string filestr,
      const std::string name_in_gmsh,
      const bool plot_to_gnuplot,
      const int step,
      const double time
      ) const;



  /*!
  \brief read restart data
  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms

  Teuchos::RCP<const Epetra_Vector> InitialGuess()  { return incvel_; }
  Teuchos::RCP<      Epetra_Vector> Residual()      { return residual_; }
  /// implement adapter fluid
  Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }
  Teuchos::RCP<const Epetra_Vector> TrueResidual()  { return trueresidual_; }
  Teuchos::RCP<const Epetra_Vector> Velnp()         { return state_.velnp_; }
  Teuchos::RCP<const Epetra_Vector> Velaf()         { return state_.velaf_; }
  Teuchos::RCP<const Epetra_Vector> Veln()          { return state_.veln_; }
  Teuchos::RCP<const Epetra_Vector> Accam()         { return state_.accam_; }
  Teuchos::RCP<const Epetra_Vector> ReadPhinp(int step);
  Teuchos::RCP<const Epetra_Vector> FsVel() {
    // get fine-scale part of velocity at time n+alpha_F or n+1
    if (Sep_ != Teuchos::null)
      AVM3Separation();

    return fsvelafStd_; } // return vector based on standard dofs

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<const DRT::DofSet> DofSet();
  Teuchos::RCP<const DRT::IndependentPBCDofSet> StandardDofSet() { return standarddofset_; }

  Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter()         { return velpressplitter_; }
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap();

  double ResidualScaling() const
  {
    if (timealgo_ == INPAR::FLUID::timeint_one_step_theta or
        timealgo_ == INPAR::FLUID::timeint_bdf2)
      return 1.0/(theta_ * dta_);
    else if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
      return alphaM_/(gamma_*dta_);
    else if (timealgo_ == INPAR::FLUID::timeint_stationary)
      return 1.0;
    else
      dserror("Not implemented so far.");
    return 0.0;
  }

  /// return time integration factor
  double TimIntParam() const;

  double EvalTime() const
  {
    return dtele_;
  };

  /// Split dof row map of discretization into velocity and pressure maps
  /*!

    We assume that each node possesses ndim velocity dofs and (optionally)
    one pressure dof. The dof row map is thus split in two.

    The velocity dofs are assigned to map 0 (the other map) and the pressure
    dofs are assigned to map 1 (the condition map).

    \author a.ger
    \date 03/08
   */
  void SetupXFluidSplit(
          const DRT::Discretization&   dis,       ///< fluid discretization
          const Teuchos::RCP<XFEM::DofManager>  dofman,    ///< information about the dofs per node
          LINALG::MapExtractor&        extractor  ///< ready made map splitter
          );

  /// Redistribute the fluid discretization and vectors according to nodegraph. This is only done in STANDARD MODE!!!
  void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph);

  /// Transfers all vectors based on an enriched dofmap to the new dofmap after a redistribution.
  void TransferVectorsToNewDistribution(Teuchos::RCP<COMBUST::FlameFront> flamefront);

  Teuchos::RCP<const Epetra_Vector> ConvectiveVel() { return Velnp(); }

  /// create field test
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();


protected:

  // don't want = operator and cctor and default constructor
  CombustFluidImplicitTimeInt operator = (const CombustFluidImplicitTimeInt& old);
  CombustFluidImplicitTimeInt (const CombustFluidImplicitTimeInt& old);
  CombustFluidImplicitTimeInt();

  /// calculate lift&drag forces and angular moments
  /*!
  Lift and drag forces are based upon the right hand side true-residual entities
  of the corresponding nodes. The contribution of the end node of a line is entirely
  added to a present L&D force.

  Idea of this routine:

  create

  std::map< label, std::set<DRT::Node*> >

  which is a set of nodes to each L&D Id
  nodal forces of all the nodes within one set are added to one L&D force

  Notice: Angular moments obtained from lift&drag forces currently refer to the
          initial configuration, i.e. are built with the coordinates X of a particular
          node irrespective of its current position.
  */
  void LiftDrag() const;

  /// helper function to get algorithm title
  std::string MapTimIntEnumToString(INPAR::FLUID::TimeIntegrationScheme term);

  /*!
  \brief give DOF information to elements
  */
  void TransferDofInformationToElements(
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      const Teuchos::RCP<XFEM::DofManager> dofmanager
      );


  /// import geometrical information about the interface and incorporate it into the fluid field
  void IncorporateInterface(const Teuchos::RCP<COMBUST::FlameFront>& flamefront);

  /*!
  \brief compute values at intermediate time steps for gen.-alpha

  */
  void GenAlphaIntermediateValues();

  /*!
  \brief update acceleration for generalized-alpha time integration

  */
  void GenAlphaUpdateAcceleration();


  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// evaluate node-based forces (currently only required for turbulent bubbly channel flow)
  void ComputeExternalForces();

  //! @name general algorithm parameters
  Teuchos::ParameterList&                       xparams_;
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle_;
  Teuchos::RCP<Epetra_Vector>                   phinp_;

  /// number of spatial dimensions
  int numdim_;

  //@}

  /// switch for type of combustion problem
  INPAR::COMBUST::CombustionType combusttype_;
  INPAR::COMBUST::VelocityJumpType veljumptype_;
  INPAR::COMBUST::FluxJumpType fluxjumptype_;

  /// do we have a turbulence model?
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  /// xfem time integration
  INPAR::COMBUST::XFEMTimeIntegration xfemtimeint_;
  INPAR::COMBUST::XFEMTimeIntegrationEnr xfemtimeint_enr_;
  INPAR::COMBUST::XFEMTimeIntegrationEnrComp xfemtimeint_enr_comp_;

  // laminar flame speed
  const double flamespeed_;
  // laminar Markstein length
  const double marksteinlength_;
  // molecular diffusivity
  const double moldiffusivity_;
  // Nitsche parameter velocity
  const double nitschevel_;
  // Nitsche parameter pressure
  const double nitschepres_;
  // boolean turning output in Gmsh files on/off
  const bool gmshoutput_;

  // parameters for two-phase flow problems
  // type of surface tension approximation
  INPAR::COMBUST::SurfaceTensionApprox surftensapprox_;
  const bool connected_interface_;

  // smoothed normal vectors for boudary integration
  const bool smoothed_boundary_integration_;
  INPAR::COMBUST::SmoothGradPhi smoothgradphi_;

  /// cpu-time measures
  double  dtele_;

  //! @name time stepping variables
  int           startsteps_; ///< maximal number of starting time steps for generalized alpha scheme
  //@}

  //! @name time step sizes
  double dtp_;    ///< time step size (used for BDF2 time integration scheme)
  //@}

  //! @name time-integration-scheme factors
  double theta_;
  double alphaM_;
  double alphaF_;
  double gamma_;
  //@}

  /// boolean indicating if stationary solution is computed as initial solution
  const bool initstatsol_;

  /// fluid reference solution iteration handling
  int itemaxFRS_; // maximal number of iterations per fgi
  int totalitnumFRS_; // total number of frs iterations in time step
  int curritnumFRS_; // number of frs iterations in fgi

  /// flag for potential nonlinear boundary conditions
  bool nonlinearbc_;

  /// pointer to turbulence statistics manager
  Teuchos::RCP<FLD::TurbulenceStatisticManager> turbstatisticsmanager_;

  /// pointer to manager for turbulent inflow transfer condition
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

  /// special discretization for xfem: holds faces and may handle extended DBC
  Teuchos::RCP<DRT::DiscretizationCombust>      xfemdiscret_;

  //! @name parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  //int dumperiod_;
  //@}

  /// flag for special flow: currently channel flow or flow in a lid-driven cavity
  std::string special_flow_;

  /// flag for fine scale subgrid viscosity used in AVM3 turbulence approach
  INPAR::FLUID::FineSubgridVisc fssgv_;

  /// flag deactivate XFEM dofs in fine-scale velocity (for AVM3)
  bool excludeXfem_;

  //! store dof distribution information for output
  Teuchos::RCP<XFEM::DofManager>      dofmanagerForOutput_;

  //! store dofset for unintersected domain
  Teuchos::RCP<DRT::IndependentPBCDofSet>                standarddofset_;

  //! store the pbc-object which maintains the PBCDofSet
  Teuchos::RCP<PeriodicBoundaryConditions> pbc_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing external loads
  Teuchos::RCP<Epetra_Vector>    external_loads_;

  /// the vector containing potential Neumann-type outflow stabilization
  Teuchos::RCP<Epetra_Vector>    outflow_stabil_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  Teuchos::RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>    incvel_;

  /// only necessary for AVM3: scale-separation matrix
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: fine-scale solution vector
  Teuchos::RCP<Epetra_Vector> fsvelafStd_;

  /// only necessary for AVM3: dof distribution map corresponding to initial unenriched configuration
  std::map<XFEM::DofKey, XFEM::DofGID>   plainnodalDofDistributionMap_;

  /// only necessary for AVM3: fine-scale solution vector
  Teuchos::RCP<Epetra_Vector> fsvelafXfem_;

  /// only necessary for AVM3: Smagorinsky constant
  double Cs_;

  /// this struct holds all state vectors needed for time integration
  ///
  /// they need to be written for restart
  struct State
  {
    /*! contains connection between nodal DofKey and a position in the global state arrays
     *
     *  should always fit to the current Epetra_Maps of this classes vectors
     */
    std::map<XFEM::DofKey, XFEM::DofGID>    nodalDofDistributionMap_;

    //! @name accelerations
    Teuchos::RCP<Epetra_Vector>    accnp_;  ///< accelerations at time n+1
    Teuchos::RCP<Epetra_Vector>    accn_;   ///< accelerations at time n
    Teuchos::RCP<Epetra_Vector>    accam_;  ///< accelerations at time n+alpha_M
    //@}

    //! @name velocities and pressures
    Teuchos::RCP<Epetra_Vector>    velnp_; ///< velocities and pressures at time n+1
    Teuchos::RCP<Epetra_Vector>    veln_;  ///< velocities and pressures at time n
    Teuchos::RCP<Epetra_Vector>    velnm_; ///< velocities and pressures at time n-1
    Teuchos::RCP<Epetra_Vector>    velaf_; ///< velocities and pressures at time n+alpha_F
    //@}
  } state_;

  /// this struct holds information about degrees of freedom for the problem
  struct PhysProb
  {
    std::set<XFEM::PHYSICS::Field> xfemfieldset_;     ///< physical XFEM fields used in this problem
    Teuchos::RCP<XFEM::ElementAnsatz> elementAnsatz_; ///< element dofs for this problem
  } physprob_;

  //! Extractor used for convergence check
  Teuchos::RCP<LINALG::MapExtractor>      velpressplitter_;

  //! Extractor used for output
  Teuchos::RCP<LINALG::MapExtractor>      velpressplitterForOutput_;

  // time integration handling
  Teuchos::RCP<XFEM::STD> timeIntStd_; //! Standard reference values
  Teuchos::RCP<XFEM::ENR> timeIntEnr_; //! Enriched reference values

  //! did we do a redistribution in this TS
  bool redist_this_step_;

  //! flag for repellant force for turbulent bubbly channel flow
  bool repellant_force_;

}; // class CombustFluidImplicitTimeInt

} // namespace FLD


#endif  // #ifndef COMBUST_FLUIDIMPLICITINTEGRATION_H
