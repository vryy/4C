/*!----------------------------------------------------------------------
\file combust_fluidimplicitintegration.H
\brief class holding implicit time integration schemes for combustion problems

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef COMBUST_FLUIDIMPLICITINTEGRATION_H
#define COMBUST_FLUIDIMPLICITINTEGRATION_H

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_io/io_ostream0.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_fluid.H"
#include "combust_flamefront.H"
#include "combust_interface.H"

#include "combust_defines.H"

namespace FLD
{
/*!
\brief implicit time integration for combustion problems

\author henke (henke@lnm.mw.tum.de)
*/
class CombustFluidImplicitTimeInt
{
  friend class CombustFluidResultTest;

public:

  /*!
  \brief constructor

  */
  explicit CombustFluidImplicitTimeInt(
      Teuchos::RCP<DRT::Discretization>  dis,
      LINALG::Solver&           solver,
      ParameterList&            params,
      IO::DiscretizationWriter& output);

  /*!
  \brief destructor

  */
  virtual ~CombustFluidImplicitTimeInt();

  /// return type of time integration scheme
  const enum INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

  /// Provide title of time integration scheme
  std::string TimIntTitle(){return MapTimIntEnumToString(TimIntScheme());}

  /// Print information about current time step to screen
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta){
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s (theta = %3f) STEP = %4d/%4d \n",
            time_,maxtime_,dta_,TimIntTitle().c_str(), theta_, step_,stepmax_);
      }
      else{
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d \n",
            time_,maxtime_,dta_,TimIntTitle().c_str(),step_,stepmax_);
      }
  }
  
  /*!
  \brief solve stationary problem
         Even though it does not make much sense, this member function belongs to the implicit time
         integration scheme for convenience. Of course, there is no time integration for stationary
         problems.
  */
  void SolveStationaryProblem();

  /// prepare a time step
  void PrepareTimeStep();

  /// prepare a nonlinear iteration
  void PrepareNonlinearSolve();

  /// import geometrical information about the interface and incorporate it into the fluid field
  void IncorporateInterface(
      const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandleNP,
      const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandleN);

  /// import the flame front from the combustion algorithm
  void StoreFlameFront(const Teuchos::RCP<COMBUST::FlameFront>& flamefront);

  /// get convection velocity vector for transfer to scalar transport field
  const Teuchos::RCP<Epetra_Vector> ConVelnp();

  /// get history vector for transfer to scalar transport field
  const Teuchos::RCP<Epetra_Vector> Hist();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or fixpoint iteration
  */
  void NonlinearSolve();

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

  /*!
  \brief lift'n'drag forces, statistics time sample and output of solution and statistics
  */
  void StatisticsAndOutput();


  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief Access output object

  */
  IO::DiscretizationWriter& DiscWriter()
  {
    return *output_;
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
      const Epetra_Map newdofrowmap);

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
  Teuchos::RCP<Epetra_Vector> CalcStresses();

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
      char* presName,
      char* velName,
      const int step,
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

  Teuchos::RCP<Epetra_Vector>& InitialGuess()  { return incvel_; }
  Teuchos::RCP<Epetra_Vector>& Residual()      { return residual_; }
  Teuchos::RCP<Epetra_Vector>& TrueResidual()  { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector>& Velnp()         { return state_.velnp_; }
  Teuchos::RCP<Epetra_Vector>& Veln()          { return state_.veln_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<DRT::Discretization>& Discretization() { return discret_; }
  Teuchos::RCP<const DRT::DofSet> StandardDofSet() { return standarddofset_; }
  Teuchos::RCP<Epetra_Vector>& RobinRHS()      { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter()         { return velpressplitter_; }
  Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.OtherMap(); }
  Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.CondMap(); }

  double Dt()   const { return dta_; }
  double Time() const { return time_; }
  int Step()    const { return step_; }
  int Itemax()  const { return itemax_; }

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  /// introduce (robin) fsi surface extractor object
  /*!
    Introduce the FSI surface to the fluid algorithm. This way special case
    algorithms (robin coupling) can work on the surface.
   */
  //void SetFSISurface(const LINALG::MapExtractor* fsisurface) { fsisurface_ = fsisurface; }

  /// calculate lift&drag forces and angular moments
  /*!
  Lift and drag forces are based upon the right hand side true-residual entities
  of the corresponding nodes. The contribution of the end node of a line is entirely
  added to a present L&D force.

  Idea of this routine:

  create

  map< label, std::set<DRT::Node*> >

  which is a set of nodes to each L&D Id
  nodal forces of all the nodes within one set are added to one L&D force

  Notice: Angular moments obtained from lift&drag forces currently refer to the
          initial configuration, i.e. are built with the coordinates X of a particular
          node irrespective of its current position.
  */
  void LiftDrag() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// integrate mass flow over a surface
  void ComputeSurfaceFlowrates() const;

  /// switch fluid field to block matrix
  void UseBlockMatrix(
      Teuchos::RCP<std::set<int> > condelements,
      const LINALG::MultiMapExtractor& domainmaps,
      const LINALG::MultiMapExtractor& rangemaps,
      bool splitmatrix=true
      );

protected:

  // don't want = operator and cctor and default constructor
  CombustFluidImplicitTimeInt operator = (const CombustFluidImplicitTimeInt& old);
  CombustFluidImplicitTimeInt (const CombustFluidImplicitTimeInt& old);
  CombustFluidImplicitTimeInt();

  /// helper function to get algorithm title
  std::string MapTimIntEnumToString(INPAR::FLUID::TimeIntegrationScheme term);

  /*!
  \brief give DOF information to elements
  */
  void TransferDofInformationToElements(
      const Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle,
      const Teuchos::RCP<XFEM::DofManager> dofmanager
      );

  //! @name general algorithm parameters
  //! fluid discretization
  Teuchos::RCP<DRT::Discretization>      discret_;
  LINALG::Solver&                        solver_;
  const ParameterList&                   params_;
  const ParameterList&                   xparams_;
  const Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! the processor ID from the communicator
  const int                     myrank_;

  //!< a singular stream to stdout in parallel environment
  IO::Ostream0 cout0_;
  //@}

  /// switch for type of combustion problem
  INPAR::COMBUST::CombustionType combusttype_;
  INPAR::COMBUST::VelocityJumpType veljumptype_;
  INPAR::COMBUST::NormalTensionJumpType normaltensionjumptype_;

  // laminar flame speed
  const double flamespeed_;
  // Nitsche parameter velocity
  const double nitschevel_;
  // Nitsche parameter pressure
  const double nitschepres_;
  // boolean indicating condensation of element dofs on/off
  const bool condensation_;

  // parameters for two-phase flow problems
  // type of surface tension approximation
  INPAR::COMBUST::SurfaceTensionApprox surftensapprox_;
  // surface tension coefficient
  const double surftenscoeff_;
  const bool connected_interface_;

  //! @name time stepping variables
  int     step_;    ///< timestep
  double  time_;    ///< physical time
  const int     stepmax_; ///< maximal number of timesteps
  const double  maxtime_; ///< maximal physical computation time
  //@}

  //! @name time step sizes
  double dta_;    ///< time step size
  double dtp_;    ///< time step size (used for BDF2 time integration scheme)
  //@}

  /// time algorithm flag
  INPAR::FLUID::TimeIntegrationScheme timealgo_;

  /// time algorithm flag for first time step
  //INPAR::FLUID::TimeIntegrationScheme startalgo_;
  
  //! @name time-integration-scheme factors
  double theta_;
//  double alphaM_;
//  double alphaF_;
//  double gamma_;
  //@}

  /// boolean indicating if stationary solution is computed as initial solution
  const bool initstatsol_;

  /// maximum number of nonlinear iterations
  int itemax_;

  /// whether we do some extrapolation at the beginning of each time step
  bool extrapolationpredictor_;

  //! @name restart variables
  const int uprestart_;
  const int upres_;
  const int writestresses_;
  //@}

  /// density extracted from element material
//  double density_;

  //! store dof distribution information for output
  Teuchos::RCP<XFEM::DofManager>      dofmanagerForOutput_;
  //! store dofset for unintersected domain
  Teuchos::RCP<DRT::DofSet>                standarddofset_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// increment of old iteration step for condensation of element dofs
  Teuchos::RCP<Epetra_Vector> oldinc_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing potential Neumann-type outflow stabilization
  Teuchos::RCP<Epetra_Vector>    outflow_stabil_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// robinrhs_: additional rhs when using robin-BC
  Teuchos::RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  Teuchos::RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>    incvel_;

  /// this struct holds all state vectors needed for time integration
  ///
  /// they need to be written for restart
  struct State
  {
    /*! contains connection between nodal DofKey and a position in the global state arrays
     *
     *  should always fit to the current Epetra_Maps of this classes vectors
     */
    std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID>    nodalDofDistributionMap_;

    /*! contains connection between elemental DofKey and a position in the global state arrays
     *
     *  should always fit to the current Epetra_Maps of this classes vectors
     */
    std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID>    elementalDofDistributionMap_;

    //! @name accelerations
    Teuchos::RCP<Epetra_Vector>    accnp_;  ///< accelerations at time n+1
    Teuchos::RCP<Epetra_Vector>    accn_;   ///< accelerations at time n
    //@}

    //! @name velocities and pressures
    Teuchos::RCP<Epetra_Vector>    velnp_;  ///< velocities and pressures at time n+1
    Teuchos::RCP<Epetra_Vector>    veln_;   ///< velocities and pressures at time n
    Teuchos::RCP<Epetra_Vector>    velnm_;  ///< velocities and pressures at time n-1
    //@}
  } state_;

  /// this struct holds information about degrees of freedom for the problem
  struct PhysProb
  {
    std::set<XFEM::PHYSICS::Field> xfemfieldset_;     ///< physical XFEM fields used in this problem
    Teuchos::RCP<XFEM::ElementAnsatz> elementAnsatz_; ///< element dofs for this problem
  } physprob_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  //! Extractor used for output
  LINALG::MapExtractor      velpressplitterForOutput_;

  //!
  Teuchos::RCP<COMBUST::FlameFront> flamefront_;

}; // class CombustFluidImplicitTimeInt

} // namespace FLD


#endif  // #ifndef COMBUST_FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
