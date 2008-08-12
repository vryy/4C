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
\brief time integration for convection-diffusion problems

\author vg (vgravem@lnm.mw.tum.de)
*/
class ScaTraTimIntImpl
{
  friend class ScaTraResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  ScaTraTimIntImpl(RCP<DRT::Discretization>       dis,
                        RCP<LINALG::Solver>            solver,
                        RCP<ParameterList>             params,
                        RCP<IO::DiscretizationWriter>  output);


  /*!
  \brief Destructor

  */
  virtual ~ScaTraTimIntImpl();


  /*!
  \brief Set the part of the righthandside belonging to the last
         timestep.
  */
  virtual void SetOldPartOfRighthandside() = 0;


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
  
  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }
  
  //! Apply Dirichlet boundary conditions on provided state vector
  void ApplyDirichletBC
  (  
    const double time,  //!< evaluation time
    Teuchos::RCP<Epetra_Vector> phinp //!< transported scalar(s)
                                      //!< (may be Teuchos::null) 
  );

  //! Apply Neumann boundary conditions
  void ApplyNeumannBC
  (
    const double time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector> phinp,   //!< transported scalar(s) state
    Teuchos::RCP<Epetra_Vector>& neumann_loads //!< Neumann loads
  );

  void Solve(bool is_stat=false);

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  virtual void Update() = 0;

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief write additional data required for restart

  */
  virtual void OutputRestart()=0;

  /*!
  \brief read restart data

  */
  virtual void ReadRestart(int step)=0;

  /*!
  \brief set the velocity field

  */
  void SetVelocityField(int veltype, int velfuncno);

  /*!
  \brief set the velocity field

  */
  void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel);

  /*!
  \brief set the initial scalar field phi

  */
  void SetInitialField(int init,int startfuncno);

  /*!
  \brief routine to calculate mass/heat flux vectors

  */
  Teuchos::RCP<Epetra_MultiVector> CalcFlux();

protected:

  // don't want = operator and cctor
  //ScaTraTimIntImpl operator = (const ScaTraTimIntImpl& old);
  //! copy constructor
  ScaTraTimIntImpl (const ScaTraTimIntImpl& old);

  RCP<DRT::Discretization>   discret_;
  RCP<LINALG::Solver>        solver_;
  RCP<ParameterList>         params_;
  RCP<IO::DiscretizationWriter>  output_;
  int                        myrank_;

  double time_;
  int step_;
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
