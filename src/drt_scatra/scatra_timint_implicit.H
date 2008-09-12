/*!----------------------------------------------------------------------
\file scatra_timint_implicit.H
\brief Associated with control routine for con.-diff.(in)stat. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_IMPLICIT_H
#define SCATRA_TIMINT_IMPLICIT_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_validparameters.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fluid/avms_solver.H"


extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../fluid_full/fluid_prototypes.h"
}

using namespace std;
using namespace Teuchos;

/// Algorithm and time integration classes for scalar transport problems
namespace SCATRA{

/*!
\brief implicit time integration for scalar transport problems

\author gjb (bauer@lnm.mw.tum.de)
*/
class ScaTraTimIntImpl
{
  friend class ScaTraResultTest;

public:

  //! Standard Constructor
  ScaTraTimIntImpl(RCP<DRT::Discretization>       dis,
                   RCP<LINALG::Solver>            solver,
                   RCP<ParameterList>             params,
                   RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntImpl();

  //! Provide enum of time integration scheme
  enum INPUTPARAMS::ScaTraTimeIntegrationScheme MethodName(){return timealgo_;}

  //! Provide title of time integration scheme
  std::string MethodTitle(){return MapTimIntEnumToString(MethodName());}

  //! Print information about current time step to screen
  void PrintTimeStepInfo()
  {
    if (myrank_==0)
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,MethodTitle().c_str(),step_,stepmax_);
    return;
  }

  //! start time loop for startingalgo, normal problems and restarts
  void Integrate();

  //! Do time integration (time loop)
  void TimeLoop();

  //! Solve stationary problem
  void SolveStationaryProblem();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  //! Apply Dirichlet boundary conditions on provided state vector
  void ApplyDirichletBC
  (
    const double& time,  //!< evaluation time
    Teuchos::RCP<Epetra_Vector>& phinp //!< transported scalar(s)
                                       //!< (may be Teuchos::null) 
  );

  //! Apply Neumann boundary conditions
  void ApplyNeumannBC
  (
    const double& time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector>& phinp,   //!< transported scalar(s) state
    Teuchos::RCP<Epetra_Vector>& neumann_loads //!< Neumann loads
  );

  void Solve(bool is_stat=false);

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  virtual void Update() = 0;

  //! update configuration and output to file
  void Output();

  //! read restart data
  virtual void ReadRestart(int step)=0;

  //! set the velocity field
  void SetVelocityField(int veltype, int velfuncno);

  //! set the velocity field
  void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel);

  //! set the initial scalar field phi
  void SetInitialField(int init,int startfuncno);

  //! set time-step-related fields for low-Mach-number flow
  void SetTimeLomaFields(RCP<const Epetra_Vector> noddensn,
                         RCP<const Epetra_Vector> noddensnm);

  //! set outer-iteration-related fields for low-Mach-number flow
  void SetIterLomaFields(int veltype,
                         RCP<const Epetra_Vector> extvel,
                         RCP<const Epetra_Vector> noddensnp);

  //! routine to return scalar field
  Teuchos::RCP<Epetra_Vector> Phinp() { return phinp_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

protected:

  // don't want = operator and cctor
  //ScaTraTimIntImpl operator = (const ScaTraTimIntImpl& old);
  //! copy constructor
  ScaTraTimIntImpl (const ScaTraTimIntImpl& old);

  //! helper function to get algorithm title
  std::string MapTimIntEnumToString
  (
    const enum INPUTPARAMS::ScaTraTimeIntegrationScheme term  //!< the enum
  );

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  //! initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep(){return;};

  //! Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  //! write state vectors (phinp and convective velocity) to BINIO
  void OutputState();

  //! write additional data required for restart
  virtual void OutputRestart() = 0;

  //! write flux vectors to BINIO
  void OutputFlux();

  //! routine to calculate mass/heat flux vectors
  Teuchos::RCP<Epetra_MultiVector> CalcFlux();

  RCP<DRT::Discretization>   discret_; /// the scalar transport discretization
  RCP<LINALG::Solver>        solver_;
  RCP<ParameterList>         params_;
  RCP<IO::DiscretizationWriter>  output_;
  int                        myrank_;

  double time_;
  int step_;
  /// flag defining type of scalar to be transported
  string scaltype_;
  int stepmax_;
  double maxtime_;

  INPUTPARAMS::ScaTraTimeIntegrationScheme timealgo_;

  int upres_;
  int uprestart_;
  string writeflux_;

  double dta_;
  double dtp_;
  double theta_;

  int cdvel_;

  /// fine-scale subgrid-viscosity flag
  string fssgd_;

  /// (standard) system matrix
  RCP<LINALG::SparseMatrix> sysmat_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>    dirichtoggle_;

  RCP<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// The residual vector
  RCP<Epetra_Vector>    residual_;

  /// phi at time n+1 and n
  RCP<Epetra_Vector>    phinp_;
  RCP<Epetra_Vector>    phin_;

  RCP<Epetra_Vector>    phinm_; // this has to go away!!

  /// density at time n+1, n and n-1
  RCP<Epetra_Vector>    densnp_;
  RCP<Epetra_Vector>    densn_;
  RCP<Epetra_Vector>    densnm_;

  /// histvector --- a linear combination of phinm, phin (BDF)
  ///                or phin, phidtn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  /// actual convective velocity
  RCP<Epetra_MultiVector>    convel_;

  /// vm3 and avms solver
  RCP<FLD::VM3_Solver> vm3_solver_;
  RCP<FLD::AVMS_Solver> avms_solver_;

  /// only necessary for vm3 solver: fine-scale subgrid-diffusivity matrix
  RCP<LINALG::SparseMatrix> sysmat_sd_;

  /// only necessary for vm3 solver: subgrid-viscosity-scaling vector
  RCP<Epetra_Vector>    sugrvisc_;


  // time measurement
  // ----------------
  double dtele_;
  double dtsolve_;

}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
#endif  // #ifdef CCADISCRET
