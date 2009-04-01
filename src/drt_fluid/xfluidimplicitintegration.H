/*!----------------------------------------------------------------------
\file xfluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef XFLUIDIMPLICITINTEGRATION_H
#define XFLUIDIMPLICITINTEGRATION_H

#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_io/io_ostream0.H"
#include "../drt_xfem/dof_management.H"

/// fluid algorithms
namespace FLD
{

/*!
\brief time integration for xfluid problems

\author (gerstenberger@lnm.mw.tum.de)
*/
class XFluidImplicitTimeInt
{
  friend class XFluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  explicit XFluidImplicitTimeInt(
      Teuchos::RCP<DRT::Discretization>  dis,
//      LINALG::Solver&                    solver,
      const ParameterList&               params
//      IO::DiscretizationWriter&          output
      );


  /*!
  \brief Destructor

  */
  virtual ~XFluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(
      const Teuchos::RCP<DRT::Discretization> cutterdiscret
      );

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(
      const Teuchos::RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
      );

  /*!
  \brief Solve stationary problem

  */
  void SolveStationaryProblem(
      const Teuchos::RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
      );

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
         Fixpoint iteration

  */
  void NonlinearSolve(
      const Teuchos::RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
      );

  /*!
  \brief (multiple) corrector

  */
  void MultiCorrector();

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.

    \param vel new guess at velocity and pressure
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);

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
  \brief lift'n'drag forces, statistics time sample and
         output of solution and statistics

  */
  void StatisticsAndOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
      Teuchos::RCP<DRT::Discretization> cutterdiscret,     ///< surface discretization used to intersect fluid elements
      int whichinitialfield,
      int startfuncno
      );

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief calculate stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<Epetra_Vector> InitialGuess() { return incvel_; }
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return state_.velnp_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return state_.veln_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter() { return velpressplitter_; }
  Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.OtherMap(); }
  Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.CondMap(); }

  double ResidualScaling() const { return density_/dta_/theta_; }
  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }
  Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

  int Itemax() const { return itemax_; }
  void SetItemax(int itemax) { itemax_ = itemax; }

  /*!
  \brief return type of time integration scheme

  */
  const FLUID_TIMEINTTYPE TimIntScheme() const {return timealgo_;}

  /// introduce free surface extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object iff we are on an ale mesh.
   */
  void SetFreeSurface(const LINALG::MapExtractor* freesurface) { freesurface_ = freesurface; }

  /// introduce (robin) fsi surface extractor object
  /*!
    Introduce the FSI surface to the fluid algorithm. This way special case
    algorithms (robin coupling) can work on the surface.
   */
  void SetFSISurface(const LINALG::MapExtractor* fsisurface) { fsisurface_ = fsisurface; }

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);
  
  /// integrate mass flow over a surface
  void ComputeSurfaceFlowRates() const;
  
  /// integrate impuls flow rate over a surface
  void ComputeSurfaceImpulsRates() const;

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool splitmatrix=true);

protected:

  // don't want = operator and cctor and default constructor
  XFluidImplicitTimeInt operator = (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt();

  /*!
  \brief compute values at intermediate time steps for gen.-alpha

  */
  void GenAlphaIntermediateValues();

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  void AssembleMatAndRHS();

  /*!
  \brief update acceleration for generalized-alpha time integration

  */
  void GenAlphaUpdateAcceleration();

  /// setup the variables to do a new nonlinear iteration
  void PrepareNonlinearSolve();

  /*!
  \brief output solution fields to Gmsh for debugging

  */
  void OutputToGmsh() const;



  /*!
  \brief compute Interface
  
    // calling this function multiple times should always result
    // in the same solution vectors and system matrix
  */
  void ComputeInterfaceAndSetDOFs(
      const Teuchos::RCP<DRT::Discretization>  cutterdiscret
      );
  
  /*!
  \brief give DOF information to elements
  */
  void TransferDofInformationToElements(
      const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
      const Teuchos::RCP<XFEM::DofManager> dofmanager
      );
  
  /*!
   * \brief plot a vector field to a Gmsh postprocessing file
   *
   * @param vectorfield      vector with solution values
   * @param filestr          basic name of the output file
   * @param name_in_gmsh     text that appears in the Gmsh window over the scale
   * @param plot_to_gnuplot  plot some values to Gnuplot (debugging - will go away sooner or later)
   */
  void PlotVectorFieldToGmsh(
      const Teuchos::RCP<Epetra_Vector>   vectorfield,
      const std::string filestr,
      const std::string name_in_gmsh,
      const bool plot_to_gnuplot
      ) const;

  //! @name general algorithm parameters
  //! fluid discretization
  const Teuchos::RCP<DRT::Discretization>      discret_;
  const Teuchos::RCP<LINALG::Solver>           solver_;
  const ParameterList&                         params_;
  const ParameterList&                         xparams_;
  const Teuchos::RCP<IO::DiscretizationWriter> output_;
  
  //! the processor ID from the communicator
  const int                           myrank_;
  //!< a singular stream to stdout in parallel environment
  IO::Ostream0 cout0_;
  
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  const bool                          alefluid_;
  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  const int     stepmax_; ///< maximal number of timesteps
  const double  maxtime_; ///< maximal physical computation time
  //@}

  /// time algorithm flag
  FLUID_TIMEINTTYPE timealgo_;
  
  /// max nonlin iter steps
  int itemax_;

  /// constant density extracted from element material for incompressible flow
  double density_;

  /// whether we do some extrapolation at the begin of each time step
  const bool extrapolationpredictor_;

  //! predictor used in nonlinear iteration
  string  predictor_;

  //! @name restart variables
  const int uprestart_;
  const int upres_;
  const int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;
  //@}

  //! @name time-integration-scheme factors
  double theta_;
  double alphaM_;
  double alphaF_;
  double gamma_;
  //@}

  //! store dof distribution information for output
  Teuchos::RCP<XFEM::DofManager>      dofmanagerForOutput_;
  //! store dofset for unintersected domain
  DRT::DofSet                dofset_out_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;
//  
//  /// (standard) system matrix
//  Teuchos::RCP<LINALG::SparseOperator> enrichsysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing potential Neumann-type outflow terms
  Teuchos::RCP<Epetra_Vector>    outflow_;

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
    std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID>   nodalDofDistributionMap_;
    
    /*! contains connection between elemental DofKey and a position in the global state arrays
     *
     *  should always fit to the current Epetra_Maps of this classes vectors
     */
    std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID>   elementalDofDistributionMap_;
    
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

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  //Teuchos::RCP<Epetra_Vector>    hist_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

//  //! Extractor used for convergence check
//  LINALG::MapExtractor      normalenrichedsplitter_;
  
  //! Extractor used for output
  LINALG::MapExtractor      velpressplitterForOutput_;

  /// free surface row dof map extractor
  const LINALG::MapExtractor* freesurface_;

  /// (robin) fsi coupling row dof map extractor
  const LINALG::MapExtractor* fsisurface_;

}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
