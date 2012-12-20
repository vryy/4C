/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o a one-step-theta time-integration scheme,

     o a two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm),

     o two variants of a generalized-alpha time-integration scheme

     and a stationary solver.

<pre>
Maintainers: Volker Gravemeier & Andreas Ehrl
             {vgravem,ehrl}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15245/-252
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H

#include "Epetra_MpiComm.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "fluid_utils_mapextractor.H" // should go away
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_fluid/fluid_timint.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"


//forward declarations

namespace DRT
{
  class Discretization;
  class ResultTest;
 // class DofSet;
 // class Condition;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace TOPOPT
{
  class Optimizer;
}
namespace ADAPTER
{
  class CouplingMortar;
}
namespace ART
{
  class ArtNetExplicitTimeInt;
}
namespace AIRWAY
{
  class RedAirwayImplicitTimeInt;
}

/// fluid algorithms
namespace FLD
{

// forward declarations
class TurbulenceStatisticManager;
class DynSmagFilter;
class LESScaleSeparation;
class Meshtying;
class TransferTurbulentInflowCondition;
class FluidMHDEvaluate;
namespace UTILS
{
  class MapExtractor;
  class KSPMapExtractor;
  class FluidImpedanceWrapper;
  class FluidWkOptimizationWrapper;
  class FluidVolumetricSurfaceFlowWrapper;
  class TotalTractionCorrector;
  class FluidInfNormScaling;
  class Fluid_couplingWrapperBase;
  template <class red_D_time_int>
  class Fluid_couplingWrapper;
}

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class FluidImplicitTimeInt : public TimInt
{
  friend class TurbulenceStatisticManager;
  friend class FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  FluidImplicitTimeInt(
      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);


  /*!
  \brief Destructor

  */
  virtual ~FluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop();

  /*!
  \brief Solve stationary problem

  */
  void SolveStationaryProblem();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
         Fixpoint iteration

  */
  void NonlinearSolve();

  /*!
  \brief solve linearised fluid

  */
  void LinearSolve();

  /*!
  \brief predictor

  */
  void Predictor();

  /*!
  \brief (multiple) corrector

  */
  void MultiCorrector();

  Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; };

  /*!
  \brief preparatives for solver

  */
  void PrepareSolve();

  /*!
  \brief update within iteration

  */
  void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment);

  /*!
 \brief convergence check

  */
  bool ConvergenceCheck(int itnum, int itmax, const double ittol);

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.
  */
  void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );

// XFEM edgebased stuff
//  void EdgeBasedStabilization(Teuchos::ParameterList & eleparams,
//      DRT::Discretization& discret,
//      int& itnum );
//
//
//
//  void findNeighborElement(DRT::ELEMENTS::Fluid * actele,
//                           RCP<DRT::Element> surface,
//                           bool & neighbor_found,
//                           int & neighbor_id);
//
//
//  void stabilizeStandardElement     ( DRT::AssembleStrategy& strategy,
//                                      DRT::ELEMENTS::Fluid *actele,
//                                      DRT::Discretization & discret);
//
//
//  void callEdgeStabandGhostPenalty(                           bool & edge_based_stab,
//                                                              bool & ghost_penalty,
//                                                              RCP<DRT::Element>       surface,
//                                                              DRT::Element * fele_1,
//                                                              DRT::Element * fele_2,
//                                                              std::vector<int> & nds_1,
//                                                              std::vector<int> & nds_2,
//                                                              DRT::Discretization &   discret,
//                                                              DRT::AssembleStrategy&  strategy);


  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

  /// Implement ADAPTER::Fluid
  void Update() { TimeUpdate(); }

  /*!
  \brief lift'n'drag forces, statistics time sample and
         output of solution and statistics

  */
  void StatisticsAndOutput();

  /*!
  \brief statistics time sample and
         output of statistics

  */
  void StatisticsOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  void OutputToGmsh(
      const int step,
      const double time,
      const bool inflow
      ) const;

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield,
    const int startfuncno
    );

  /// Implement ADAPTER::Fluid
  Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

  /// Implement ADAPTER::Fluid
  Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(Teuchos::RCP<const Epetra_Vector> velpres);

  /*!
  \brief reset density in topology optimization problems

  */
  void SetTopOptData(
      RCP<const Epetra_Vector> porosity,
      RCP<TOPOPT::Optimizer>& optimizer
  );

  /// Reset state vectors
  void Reset(
      bool completeReset = false,
      bool newFiles = false,
      int iter = -1);

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief evaluate divergence of velocity field

  */
  void EvaluateDivU();

  /*!
  \brief calculate stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief calculate wall shear stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcWallShearStresses();

  /*!
  \brief calculate subfilter stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcSFS(
  const int   i,
  const int   j);

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);
  /*!
  \brief get restart data in case of turbulent inflow computation

  */
  void SetRestart(const int step,
                  const double time,
                  Teuchos::RCP<const Epetra_Vector> readvelnp,
                  Teuchos::RCP<const Epetra_Vector> readveln,
                  Teuchos::RCP<const Epetra_Vector> readvelnm,
                  Teuchos::RCP<const Epetra_Vector> readaccnp,
                  Teuchos::RCP<const Epetra_Vector> readaccn);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<const Epetra_Vector> InitialGuess() { return incvel_; }

  /// return implemented residual (is not an actual force in Newton [N])
  Teuchos::RCP<      Epetra_Vector> Residual() { return residual_; }

  /// implement adapter fluid
  Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }

  /// Return true residual, ie the actual force in Newton [N]
  Teuchos::RCP<const Epetra_Vector> TrueResidual() { return trueresidual_; }

  Teuchos::RCP<const Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<      Epetra_Vector> ViewOfVelnp() { return velnp_; }
  Teuchos::RCP<const Epetra_Vector> Velaf() { return velaf_; }
  Teuchos::RCP<const Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<const Epetra_Vector> Velnm() { return velnm_; }
  Teuchos::RCP<const Epetra_Vector> Accnp() { return accnp_; }
  Teuchos::RCP<const Epetra_Vector> Accn() { return accn_; }
  Teuchos::RCP<const Epetra_Vector> Accam() { return accam_; }
  Teuchos::RCP<const Epetra_Vector> Scaaf() { return scaaf_; }
  Teuchos::RCP<const Epetra_Vector> Scaam() { return scaam_; }
  Teuchos::RCP<const Epetra_Vector> Hist() { return hist_; }
  Teuchos::RCP<const Epetra_Vector> GridVel() { return gridv_; }
  Teuchos::RCP<      Epetra_Vector> ViewOfGridVel() { return gridv_; }
  Teuchos::RCP<const Epetra_Vector> FsVel() { return fsvelaf_; }

  /// provide access to the Dirichlet map
  Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// implement ADAPTER::Fluid, pass through DirichMaps()
  Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return DirichMaps(); }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Shorten the Dirichlet DOF set
  ///
  /// The method shortens the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the monolithic FSI where fluid boundary conditions on the
  /// FSI interface have to be removed if it is a fluid split.
  virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  Teuchos::RCP<const Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<      Epetra_Vector> ViewOfDispnp() { return dispnp_; }
  Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter(){ return velpressplitter_; };
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap();
//  void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);
//  double TimeScaling() const;

  /// Use ResidualScaling() to convert the implemented fluid residual to an actual force with unit Newton [N]
  /*! In order to avoid division by time step size \f$\Delta t\f$
   *  the fluid balance of linear momentum is implemented in a way
   *  that the residual does not have the unit Newton [N].
   *  By multiplication with ResidualScaling() the residual is
   *  converted to the true residual in unit Newton [N], ie a real force.
   *
   *  \sa trueresidual_
   *  \sa TrueResidual()
   */
  double ResidualScaling() const
  {
    if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
      return 1.0;
    else if (is_genalpha_)
      return alphaM_/(gamma_*dta_);
    else
      return 1.0/(theta_*dta_);
  }

  /// return time integration factor
  double TimIntParam() const;

  /// update velocity increment after Newton step
  virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel);

//  int Itemax() const { return params_->get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

  void PrepareTimeandStep() { UpdateStepandTime(); }
  void UpdateStepandTime(){ step_ += 1;
  time_ += dta_;}


  /*!
  \brief return type of time integration scheme

  */
//  const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

  /*!
  \brief set fields for low-Mach-number flow within outer iteration loop

  */
  void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                         RCP<const Epetra_Vector> scalaram,
                         RCP<const Epetra_Vector> scalardtam,
                         RCP<const Epetra_Vector> fsscalaraf,
                         const double             thermpressaf,
                         const double             thermpressam,
                         const double             thermpressdtaf,
                         const double             thermpressdtam,
                         Teuchos::RCP<DRT::Discretization> scatradis);

  /*!
  \brief set fields for low-Mach-number flow at end of time step

  */
  void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                         const double             thermpressnp,
                         RCP<const Epetra_Vector> scatraresidual,
                         Teuchos::RCP<DRT::Discretization> scatradis,
                         const int                whichscalar = -1);

  /*!
  \brief set velocity field obtained by separate computation

  */
  void SetVelocityField(Teuchos::RCP<const Epetra_Vector> setvelnp) { velnp_->Update(1.0,*setvelnp,0.0); return;}

  /// provide access to turbulence statistics manager
  Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager();
  /// provide access to the box filter for dynamic Smagorinsky model
  Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter();

  /// introduce surface split extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object if we are on an ale mesh.
   */
  void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) { surfacesplitter_ = surfacesplitter; }

  /// determine grid velocity
  void UpdateGridv();

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;

  /// compute flow rate
  void ComputeFlowRates() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool splitmatrix=true);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

  //! @name methods for turbulence models

  void ApplyScaleSeparationForLES();

  void OutputofFilteredVel(Teuchos::RCP<Epetra_Vector> outvec,
                           Teuchos::RCP<Epetra_Vector> fsoutvec);

  void PrintTurbulenceModel();
  //@}

  /// set the initial porosity field
  void SetInitialPorosityField(
    const INPAR::POROELAST::InitialField,   ///< type of initial field
	//const int,                                ///< type of initial field
    const int startfuncno                   ///< number of spatial function
    );

  void UpdateIterIncrementally(
    Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
    );

  //! @name methods for fsi
  /// Extrapolation of vectors from mid-point to end-point t_{n+1}
  Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
      Teuchos::RCP<Epetra_Vector> vecn, ///< vector at time level t_n
      Teuchos::RCP<Epetra_Vector> vecm  ///< vector at time level of equilibrium
      );
  //@}

  /// create field test
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  Teuchos::RCP<const Epetra_Vector> ConvectiveVel();

  /*! \brief Calculate a integrated divergence operator in vector form
   *
   *   The vector valued operator \f$B\f$ is constructed such that
   *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
   */
  Teuchos::RCP<Epetra_Vector> CalcDivOp();

protected:

  // don't want = operator and cctor
  FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

  /*!
  \brief  increment time and step value

  */
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  void AssembleMatAndRHS();

  /*!
  \brief apply Dirichlet boundary conditions to system of equations

  */
  void ApplyDirichletToSystem();

  /*!
  \brief update acceleration for generalized-alpha time integration

  */
  void GenAlphaUpdateAcceleration();

  /*!
  \brief compute values at intermediate time steps for gen.-alpha

  */
  void GenAlphaIntermediateValues();


  //! @name general algorithm parameters
  //! project nullspace (useful for completely Dirichlet bounded problems)
  bool                          project_;
  //! integrated basis vector (pressure component)
  Teuchos::RCP<Epetra_Vector>   w_;
  //! basis vector for nullspace (vector of ones in the pressure component)
  Teuchos::RCP<Epetra_Vector>   c_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  //! do we have a turbulence model?
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  //@}

  /// number of spatial dimensions
  int numdim_;

  //! @name time stepping variables
  int     numstasteps_; ///< number of steps for starting algorithm
  bool    startalgo_;   ///< flag for starting algorithm
  //@}

  // genalpha integration scheme (afgenalpha or npgenalpha)
  bool is_genalpha_;

  /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
  /// and at n+alpha_M/n as well as its time derivative at n+alpha_F/n+1 and n+alpha_M/n
  double thermpressaf_;
  double thermpressam_;
  double thermpressdtaf_;
  double thermpressdtam_;

  /// gas constant (only for low-Mach-number flow)
  double gasconstant_;

  //! use (or not) linearisation of reactive terms on the element
  INPAR::FLUID::LinearisationAction newton_;

  /// whether we do some extrapolation at the begin of each time step (additional switch)
  bool extrapolationpredictor_;

  //! kind of predictor used in nonlinear iteration
  string  predictor_;

  //! @name restart variables
  int writestresses_;
  int write_wall_shear_stresses_;
  //@}

  //! @name time step sizes
  double dtp_;
  //@}

  //! @name time-integration-scheme factors
  double theta_;
  double omtheta_;
  double alphaM_;
  double alphaF_;
  double gamma_;
  //@}

  //! @name parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;
  //@}

  //! @name norms for convergence check
  double incvelnorm_L2_;
  double incprenorm_L2_;
  double velnorm_L2_;
  double prenorm_L2_;
  double vresnorm_;
  double presnorm_;
  //@}

  /// flag for special flow: currently channel flow or flow in a lid-driven cavity
  string special_flow_;

  /// flag for potential Neumann inflow
  bool neumanninflow_;

  /// flag for potential 3D Reduced_D coupling
  bool strong_redD_3d_coupling_;

  /// form of convective term
  string convform_;

  /// fine-scale subgrid-viscosity flag
  string fssgv_;

  /// cpu-time measures
  double  dtele_;
  double  dtfilter_;
  double  dtsolve_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// linearization with respect to mesh motion
  RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<FLD::UTILS::MapExtractor> vol_flow_rates_bc_extractor_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing potential Neumann-type outflow terms
//  RCP<Epetra_Vector>    outflow_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;
  RCP<Epetra_Vector>    trac_residual_;

  /// rhs: right hand side vector for linearised flow solution
//  RCP<Epetra_Vector>    rhs_;

  /// robinrhs_: additional rhs when using robin-BC
  RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions (Unit: Newton [N])
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
  RCP<Epetra_Vector>    accnp_;
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accam_;
  //@}

  //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  RCP<Epetra_Vector>    velaf_;
  //@}

  //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
  RCP<Epetra_Vector>    scaaf_;
  RCP<Epetra_Vector>    scaam_;
  //@}

  //! @name displacements at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;
  //@}

  //! @name optimization data

  /// density in topology optimization
  RCP<const Epetra_Vector>    topopt_porosity_;

  /// optimizer for topology optimization
  RCP<TOPOPT::Optimizer> optimizer_;
  //@}

  /// only necessary for AVM3: scale-separation matrix
  RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: fine-scale solution vector
  RCP<Epetra_Vector>    fsvelaf_;

  /// only necessary for Scale Similarity Model: filtered velocity
  RCP<Epetra_MultiVector>    filteredvel_;

  /// only necessary for Scale Similarity Model: filtered reynolds stress
  RCP<Epetra_MultiVector>    filteredreystr_;

  /// only necessary for Scale Similarity Model: fine scale velocity
  RCP<Epetra_MultiVector>    finescalevel_;

  /// only necessary for LES models including filtered quantities: filter type
  enum INPAR::FLUID::ScaleSeparation scale_sep_;

  /// fine-scale scalar: only necessary for multifractal subgrid-scale modeling in loma
  RCP<Epetra_Vector>    fsscaaf_;

  /// grid velocity (set from the adapter!)
  RCP<Epetra_Vector>    gridv_;

  /// initial porosity (poroelasticity)
  RCP<Epetra_Vector>    initporosityfield_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;


  //! manager for turbulence statistics
  RCP<FLD::TurbulenceStatisticManager> statisticsmanager_;

  //! @name Dynamic Smagorinsky model: methods and variables
  //        -------------------------

  //! one instance of the filter object
  RCP<FLD::DynSmagFilter> DynSmag_;
  RCP<LESScaleSeparation> ScaleSepGMO_;

  //@}

  //! @name Mixed-hybrid Dirichlet condition
  //        --------------------------------
  //
  //! object for a redistributed evaluation of the mixed-hybrid Dirichlet condition
  RCP<FLD::FluidMHDEvaluate> MHD_evaluator_;
  //@}

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// ALE dof map
  LINALG::MapExtractor                   meshmap_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;

  /// Krylov space projection map extractor
  FLD::UTILS::KSPMapExtractor kspsplitter_;

  /// connection between master and slave nodes on this proc (based on row map!!!!)
  RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave_;

  /// a manager doing the transfer of boundary data for
  /// turbulent inflow profiles from a separate (periodic) domain
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

  /// bio related special (outflow) boundaries
  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  /// bio related special (outflow) optimization of boundaries
  RCP<UTILS::FluidWkOptimizationWrapper> Wk_optimization_;

  /// bio related special (in/outflow) boundaries
  RCP<UTILS::FluidVolumetricSurfaceFlowWrapper> vol_surf_flow_bc_;

  /// bio related special (in/outflow) traction velocity component adder
  RCP<UTILS::TotalTractionCorrector> traction_vel_comp_adder_bc_;

  /// bio related, 3D to reduced-D coupling
  RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_art_;

  /// 1D arterial network time integration
  RCP<ART::ArtNetExplicitTimeInt>   ART_exp_timeInt_;

  /// bio related, 3D to reduced-D coupling
  RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_airways_;

  /// 1D arterial network time integration
  RCP<AIRWAY::RedAirwayImplicitTimeInt>   airway_imp_timeInt_;

  /// @name special relaxation state

  bool inrelaxation_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  /// maps for womersley flow profile which is applied as a Dirichlet condition
  RCP<Epetra_Map>                    vol_surf_flow_bcmaps_;

  /// coupling of fluid-fluid at an internal interface
  RCP<FLD::Meshtying>   meshtying_;

  /// flag for mesh-tying
  int msht_;

  //@}

  // possible inf-norm scaling of linear system / fluid matrix
  RCP<FLD::UTILS::FluidInfNormScaling> fluid_infnormscaling_;

private:

  void SplitMatrixAndVector(RCP<Epetra_Vector>    incvel_);

  //@}

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementGeneralFluidParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTimeParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTurbulenceParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementLomaParameter();

  //! @name Set general parameter in class f3Parameter for porous media problems
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementPoroParameter();

  //! @name Special method for turbulent variable-density flow at low Mach number with multifractal subgrid-scale modeling
  /*!

  \brief adaptation of CsgsD to CsgsB
         Since CsgsB depends on the resolution if the near-wall limit is included,
         CsgsD is adapted accordingly by using the mean value of the near-wall  correction.

  */
  void RecomputeMeanCsgsB();

}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
