/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

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

using namespace std;
using namespace Teuchos;

//forward declarations

namespace DRT
{
  class Discretization;
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

namespace ADAPTER
{
  class CouplingMortar;
}
#ifdef D_ARTNET
namespace ART
{
  class ArtNetExplicitTimeInt;
}
#endif //D_ARTNET

#ifdef D_RED_AIRWAYS
namespace AIRWAY
{
  class RedAirwayImplicitTimeInt;
}
#endif //D_RED_AIRWAYS

/// fluid algorithms
namespace FLD
{

// forward declarations
class TurbulenceStatisticManager;
class DynSmagFilter;
class LESScaleSeparation;
class Meshtying;
class TransferTurbulentInflowCondition;
namespace UTILS
{
  class MapExtractor;
  class KSPMapExtractor;
  class FluidImpedanceWrapper;
  class FluidWkOptimizationWrapper;
  class FluidVolumetricSurfaceFlowWrapper;
  class TotalTractionCorrector;
#ifdef D_ARTNET
  class Fluid_couplingWrapperBase;
  template <class red_D_time_int>
  class Fluid_couplingWrapper;
#else
#ifdef D_RED_AIRWAYS
  class Fluid_couplingWrapperBase;
  template <class red_D_time_int>
  class Fluid_couplingWrapper;
#endif
#endif
}

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class FluidImplicitTimeInt
{
  friend class TurbulenceStatisticManager;
  friend class FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  FluidImplicitTimeInt(RCP<DRT::Discretization>  actdis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


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

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

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
  \brief Access output object

  */
  IO::DiscretizationWriter& DiscWriter()
  {
    return output_;
  }

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield,
    const int startfuncno
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
   * \brief calculate fluid volume (2d case: calculate fluid "surface")
   */
  double CalcFluidVolume();

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

  Teuchos::RCP<Epetra_Vector> InitialGuess() { return incvel_; }
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Velaf() { return velaf_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> Velnm() { return velnm_; }
  Teuchos::RCP<Epetra_Vector> Accnp() { return accnp_; }
  Teuchos::RCP<Epetra_Vector> Accn() { return accn_; }
  Teuchos::RCP<Epetra_Vector> Accam() { return accam_; }
  Teuchos::RCP<Epetra_Vector> Scaaf() { return scaaf_; }
  Teuchos::RCP<Epetra_Vector> Scaam() { return scaam_; }
  Teuchos::RCP<Epetra_Vector> Hist() { return hist_; }
  Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<Epetra_Vector> Dispn() { return dispn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter(){ return velpressplitter_; };
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap();

  double ResidualScaling() const
  { if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
         return 1.0;
    else if (TimIntScheme()==INPAR::FLUID::timeint_afgenalpha)
         return alphaM_/(gamma_*dta_);
    else return 1.0/(theta_*dta_);}
  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  void UpdateStepandTime(){ step_ += 1;
  time_ += dta_;}


  /*!
  \brief return type of time integration scheme

  */
  const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

  /*!
  \brief set fields for low-Mach-number flow within outer iteration loop

  */
  void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                         RCP<const Epetra_Vector> scalaram,
                         RCP<const Epetra_Vector> scalardtam,
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

  /// introduce surface split extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object iff we are on an ale mesh.
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


protected:

  // don't want = operator and cctor
  FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

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
  //! fluid discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  //! project nullspace (useful for completely Dirichlet bounded problems)
  bool                          project_;
  //! integrated basis vector (pressure component)
  Teuchos::RCP<Epetra_Vector>   w_;
  //! basis vector for nullspace (vector of ones in the pressure component)
  Teuchos::RCP<Epetra_Vector>   c_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  //! do we have a turblence model?
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  //@}

  /// number of spatial dimensions
  int numdim_;

  //! @name time stepping variables
  double  time_;        ///< physical time
  int     step_;        ///< timestep
  int     stepmax_;     ///< maximal number of timesteps
  double  maxtime_;     ///< maximal physical computation time
  int     numstasteps_; ///< number of steps for starting algorithm
  bool    startalgo_;   ///< flag for starting algorithm
  //@}

  /// time algorithm flag
  INPAR::FLUID::TimeIntegrationScheme timealgo_;
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
  int uprestart_;
  int upres_;
  int writestresses_;
  int write_wall_shear_stresses_;
  //@}

  //! @name time step sizes
  double dta_;
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

  /// flag for physical type of fluid flow (standard: incompressible)
  enum INPAR::FLUID::PhysicalType physicaltype_;

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
  RCP<Epetra_Vector>    outflow_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;
  RCP<Epetra_Vector>    trac_residual_;

  /// rhs: right hand side vector for linearised flow solution
  RCP<Epetra_Vector>    rhs_;

  /// robinrhs_: additional rhs when using robin-BC
  RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
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

  //! @name displacemets at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;
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

  /// grid velocity (set from the adapter!)
  RCP<Epetra_Vector>    gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;


  //! manager for turbulence statistics
  RCP<FLD::TurbulenceStatisticManager> statisticsmanager_;

  //! @name Dynamic Smagorinsky model: methods and variables
  //        -------------------------

  //! one instance of the filter object
  RCP<DynSmagFilter> DynSmag_;
  RCP<LESScaleSeparation> ScaleSepGMO_;

  //@}

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;

  /// Krylov space projection map extractor
  FLD::UTILS::KSPMapExtractor kspsplitter_;

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > > pbcmapmastertoslave_;

  /// a manager doing the transfer of boundary data for
  /// turbulent inflow profiles from a separate (periodic) domain
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

  /// bio related special (outflow) boudaries
  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  /// bio related special (outflow) optimization of boundaries
  RCP<UTILS::FluidWkOptimizationWrapper> Wk_optimization_;

  /// bio related special (in/outflow) boudaries
  RCP<UTILS::FluidVolumetricSurfaceFlowWrapper> vol_surf_flow_bc_;

  /// bio related special (in/outflow) traction velocity component adder
  RCP<UTILS::TotalTractionCorrector> traction_vel_comp_adder_bc_;

#ifdef D_ARTNET
  /// bio related, 3D to redued-D coupling
  RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_art_;

  /// 1D arterial network time integration
  RCP<ART::ArtNetExplicitTimeInt>   ART_exp_timeInt_;
#endif // D_ARTNET

#ifdef D_RED_AIRWAYS

  /// bio related, 3D to redued-D coupling
  RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_airways_;

  /// 1D arterial network time integration
  RCP<AIRWAY::RedAirwayImplicitTimeInt>   airway_imp_timeInt_;

#endif //D_RED_AIRWAYS
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

}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
