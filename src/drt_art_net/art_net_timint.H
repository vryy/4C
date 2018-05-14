/*!----------------------------------------------------------------------
\file art_net_timint.H
\brief Control routine for arterial network (time) integration.

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/

#ifndef SRC_DRT_ART_NET_ART_NET_TIMINT_H_
#define SRC_DRT_ART_NET_ART_NET_TIMINT_H_

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_art_net.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_bio.H"


/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace LINALG
{
  class Solver;
}

namespace DRT
{
  class ResultTest;
}

/*==========================================================================*/
//! Algorithm and time integration classes for artery network problems
/*==========================================================================*/
namespace ART
{

/*!
 * \brief time integration for artery network problems
 */

class TimInt : public ADAPTER::ArtNet
{

public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  //! Standard Constructor
  TimInt(Teuchos::RCP<DRT::Discretization>  dis,
         const int linsolvernumber,
         const Teuchos::ParameterList&   probparams,
         const Teuchos::ParameterList&   artparams,
         FILE* errfile,
         IO::DiscretizationWriter&  output);

  //! Destructor
  virtual ~TimInt();

  //! initialize time integration
  virtual void Init(
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& arteryparams,
      const   std::string& scatra_disname
      );

  //! get discretization
  Teuchos::RCP<DRT::Discretization > Discretization() { return discret_; }

  virtual double  Dt()   const { return dta_; }

  double  Time() const { return time_; }
  int     Step() const { return step_; }

  int     Itemax() const { return params_.get<int>("max nonlin iter steps"); }

  virtual void Output(bool               CoupledTo3D,
      Teuchos::RCP<Teuchos::ParameterList> CouplingParams) = 0;

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

  /*!
  \brief prepare the loop

  */
  virtual void PrepareTimeLoop();

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

  //! set the initial field on the artery discretization
  virtual void SetInitialField(
    const INPAR::ARTDYN::InitialField  init,       //!< type of initial field
    const int                          startfuncno //!< number of spatial function
  )
  {
    // each artery integration should overwrite this if used
    dserror("not implemented");
  }


  /// setup the variables to do a new time step
  virtual void PrepareTimeStep();

  /// setup the variables to do a new time step
  virtual void PrepareLinearSolve()
  {
    // each artery integration should overwrite this if used
    dserror("not implemented");
  }

  /// setup the variables to do a new time step
  virtual void AssembleMatAndRHS()
  {
    // each artery integration should overwrite this if used
    dserror("not implemented");
  }

  /// direct access to system matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
  {
    return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
  };

  //! right-hand side alias the dynamic force residual
  virtual Teuchos::RCP<const Epetra_Vector> RHS() const
  {
    return rhs_;
  }

  //! iterative update of primary variable
  virtual void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
  {
    // each artery integration should overwrite this if used
    dserror("not implemented");
    return;
  }

  // get solution vector
  virtual Teuchos::RCP<const Epetra_Vector> Pressurenp() const
  {
    // each artery integration should overwrite this if used
    dserror("not implemented");
    return Teuchos::null;
  }
  /*!
  \brief solve linearised artery and bifurcation

  */
  virtual void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) = 0;

  virtual void SolveScatra() = 0;

  //! is output needed for the current time step?
  bool DoOutput()
    { return ((step_%upres_==0 )or (step_%uprestart_==0)); };

  // set solve scatra flag
  void SetSolveScatra(const bool solvescatra){ solvescatra_ = solvescatra; return;}

  //! Return MapExtractor for Dirichlet boundary conditions
  virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const
  {
   return dbcmaps_;
  }

  // create field test
  virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;

protected:

  //! @name general algorithm parameters
  //! arterial network discretization
  Teuchos::RCP<DRT::Discretization>      discret_;
  //! linear solver
  Teuchos::RCP<LINALG::Solver>         solver_;
  const Teuchos::ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! error file handle
  FILE* errfile_;
  //! the processor ID from the communicator
  int                           myrank_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// rhs: right hand side vector
  Teuchos::RCP<Epetra_Vector>    rhs_;

  /// (scatra) system matrix
  Teuchos::RCP<LINALG::SparseOperator> scatra_sysmat_;

  /// rhs: right hand side vector of scatra
  Teuchos::RCP<Epetra_Vector>    scatra_rhs_;

  //! @name time step sizes
  double dta_;
  double dtp_;

  //! @name cpu-time measures
  double  dtele_;
  double  dtsolve_;
  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}

  //! @name restart variables
  int uprestart_;
  int upres_;
  //@}

  bool solvescatra_;
  const int linsolvernumber_;

  //!
  bool coupledTo3D_;
  //@}


}; // class TimInt
} // namespace ArtNet




#endif /* SRC_DRT_ART_NET_ART_NET_TIMINT_H_ */
