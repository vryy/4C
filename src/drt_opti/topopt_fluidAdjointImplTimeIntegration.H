/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointImplTimeIntegration.H

\brief fluid adjoint implicit time integration for topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_
#define TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_


#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_topopt.H"


/// forward declarations
class Epetra_Map;
class Epetra_Vector;

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
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace TOPOPT
{
  class Optimizer;
}

namespace TOPOPT
{
namespace ADJOINT
{
/*
  \brief Associated with control routine for adjoint fluid (in)stationary solvers,

     including instationary solvers based on
     - one-step-theta time-integration scheme

     and stationary solver.

    \author winklmaier
    \date 02/12
 */
class ImplicitTimeInt
{
//  friend class FLD::FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  ImplicitTimeInt(Teuchos::RCP<DRT::Discretization>  actdis,
                       LINALG::Solver&           solver,
                       Teuchos::ParameterList&            params,
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
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief update configuration and output to gmsh file

  */  void OutputToGmsh(
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
  void SetInitialAdjointField(
    const INPAR::TOPOPT::InitialAdjointField initfield,
    const int startfuncno
    );

  /*!
  \brief reset density in topology optimization problems

  */
  void SetTopOptData(
      Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluidvelocities,
      Teuchos::RCP<Epetra_Vector> porosity,
      Teuchos::RCP<TOPOPT::Optimizer> optimizer
  );

  void EvaluateDirichlet();

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter(){ return velpressplitter_; };
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap();

  double ResidualScaling() const;

  double Dt() const { return dt_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }


  /*!
  \brief return type of time integration scheme

  */
  const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

  /*!
  \brief set velocity field obtained by separate computation

  */
  void SetAdjointField(Teuchos::RCP<const Epetra_Vector> setvelnp);

protected:

  // don't want = operator and cctor
  ImplicitTimeInt operator = (const ImplicitTimeInt& old);
  ImplicitTimeInt (const ImplicitTimeInt& old);

  /*!
  \brief  increment time and step value

  */
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ -= dt_;
  }


  //! @name general algorithm parameters
  //! fluid discretization
  const Teuchos::RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  Teuchos::ParameterList&       params_;
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
  //@}

  /// time algorithm flag
  INPAR::FLUID::TimeIntegrationScheme timealgo_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  //@}

  //! @name time step size
  double dt_;
  //@}

  //! @name time-integration-scheme factors
  double theta_;
  double omtheta_;
  double theta_pre_;
  double omtheta_pre_;
  double theta_div_;
  double omtheta_div_;
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
  std::string special_flow_;

  /// flag for potential Neumann inflow
  bool neumanninflow_;

  /// cpu-time measures
  double  dtele_;
  double  dtsolve_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>    incvel_;

  //! @name velocity and pressure at time n+1, n, n-1
  Teuchos::RCP<Epetra_Vector>    velnp_;
  Teuchos::RCP<Epetra_Vector>    veln_;
  //@}

  /// map <timestep,fluid velocity/pressure>
  Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluid_vels_;

  /// density in topology optimization
  Teuchos::RCP<Epetra_Vector>    topopt_porosity_;

  /// optimizer for topology optimization
  Teuchos::RCP<TOPOPT::Optimizer> optimizer_;

  //! Extractor used for convergence check
  Teuchos::RCP<LINALG::MapExtractor>      velpressplitter_;

  /// connection between master and slave nodes on this proc
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave_;

private:

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
