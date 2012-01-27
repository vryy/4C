/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointImplTimeIntegration.H

\brief Associated with control routine for adjoint fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_
#define TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_


#include <Epetra_Vector.h>
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H" // should go away

using namespace std;
using namespace Teuchos;

//forward declarations
namespace DRT
{
  class Discretization;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
//  class MultiMapExtractor;
//  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace TOPOPT
{
  class Optimizer;
}

/// fluid algorithms
namespace FLD
{
  namespace UTILS
  {
    class MapExtractor;
    class KSPMapExtractor;
    class FluidImpedanceWrapper;
    class FluidWkOptimizationWrapper;
    class FluidVolumetricSurfaceFlowWrapper;
    class TotalTractionCorrector;
  }
}

namespace TOPOPT
{
namespace ADJOINT
{
/*!
\brief time integration for adjoint fluid problems in topology optimization

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class ImplicitTimeInt
{
//  friend class FLD::FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  ImplicitTimeInt(RCP<DRT::Discretization>  actdis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~ImplicitTimeInt();

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
  void NonLinearSolve();

  /*!
  \brief predictor

  */
  void Predictor();

  /*!
  \brief preparatives for solver

  */
  void PrepareSolve();

  /*!
  \brief update within iteration

  */
  void IterUpdate();

  /*!
 \brief convergence check

  */
  void ConvergenceCheck();

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
  \brief reset density in topology optimization problems

  */
  void SetTopOptData(
      RCP<Epetra_Vector> porosity,
      RCP<TOPOPT::Optimizer> optimizer
  );

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

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

  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Velaf() { return velaf_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> Velnm() { return velnm_; }
  Teuchos::RCP<Epetra_Vector> Accnp() { return accnp_; }
  Teuchos::RCP<Epetra_Vector> Accn() { return accn_; }
  Teuchos::RCP<Epetra_Vector> Accam() { return accam_; }
  Teuchos::RCP<Epetra_Vector> Hist() { return hist_; }

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
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

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
  \brief set velocity field obtained by separate computation

  */
  void SetVelocityField(Teuchos::RCP<const Epetra_Vector> setvelnp) { velnp_->Update(1.0,*setvelnp,0.0); return;}

protected:

  // don't want = operator and cctor
  ImplicitTimeInt operator = (const ImplicitTimeInt& old);
  ImplicitTimeInt (const ImplicitTimeInt& old);

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
  const RCP<DRT::Discretization>      discret_;
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

  //! use (or not) linearisation of reactive terms on the element
  INPAR::FLUID::LinearisationAction newton_;

  /// whether we do some extrapolation at the begin of each time step (additional switch)
  bool extrapolationpredictor_;

  //! @name restart variables
  int uprestart_;
  int upres_;
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

  /// flag for special flow: currently channel flow or flow in a lid-driven cavity
  string special_flow_;

  /// flag for potential Neumann inflow
  bool neumanninflow_;

  /// form of convective term
  string convform_;

  /// cpu-time measures
  double  dtele_;
  double  dtsolve_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

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

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

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

  //! @name displacemets at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;
  //@}

  /// density in topology optimization
  RCP<Epetra_Vector>    topopt_porosity_;

  /// optimizer for topology optimization
  RCP<TOPOPT::Optimizer> optimizer_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// row dof map extractor
  const FLD::UTILS::MapExtractor* surfacesplitter_;

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > > pbcmapmastertoslave_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

private:

  void SplitMatrixAndVector(RCP<Epetra_Vector>    incvel_);

  //@}

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementGeneralAdjointParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTimeParameter();

}; // class FluidImplicitTimeInt
}
}



#endif /* TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_ */
#endif  // #ifdef CCADISCRET
