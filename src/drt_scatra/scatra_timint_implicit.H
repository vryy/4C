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
#include "../drt_inpar/drt_validparameters.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fluid/avms_solver.H"
#include "../drt_fluid/drt_periodicbc.H"

#include "../drt_inpar/inpar_scatra.H"


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
  enum INPAR::SCATRA::TimeIntegrationScheme MethodName(){return timealgo_;}

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

  //! do time integration (time loop)
  void TimeLoop(bool nonlinear = false);

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  //! contains the nonlinear iteration loop
  void NonlinearSolve();

  //! check convergence (or divergence) of nonlinear iteration
  bool AbortNonlinIter(const int itnum, const int itmax, const double ittol, double& actresidual);

  //! solve linear system
  void Solve();

  //! update the solution after convergence of the nonlinear iteration.
  virtual void Update() = 0;

  //! update configuration and output to file
  void Output();

  //! read restart data
  virtual void ReadRestart(int step)=0;

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  /// note: VM3 solver still needs an explicit toggle vector for construction
  const Teuchos::RCP<const Epetra_Vector> DirichletToggle();

  //! prepare AVM3-based scale separation
  void AVM3Preparation();

  //! scaling of AVM3-based subgrid-diffusivity matrix
  void AVM3Scaling(ParameterList& eleparams);

  //! set the velocity field
  void SetVelocityField(int veltype, int velfuncno);

  //! set the velocity field
  void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel);

  //! set the initial scalar field phi
  void SetInitialField(int init,int startfuncno);

  //! set velocity field for low-Mach-number flow
  virtual void SetLomaVelocity(RCP<const Epetra_Vector> extvel,
                       RCP<DRT::Discretization> fluiddis)=0;

  //! predict density field for low-Mach-number flow
  virtual void PredictDensity() = 0;

  //! compute density field for low-Mach-number flow
  void ComputeDensity();

  //! compute time derivative of density field for low-Mach-number flow
  virtual void ComputeDensityDerivative() = 0;

  //! check convergence of density field for low-Mach-number flow
  bool DensityConvergenceCheck(int itnum, int itmax, const double ittol);

  //! update density field for low-Mach-number flow
  virtual void UpdateDensity() = 0;

  //! output of mean temperature and density
  void OutputMeanTempAndDens();

  //! routine to return density field at time step n+1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> DensNp() { return densnp_; }

  //! routine to return density field at time step n for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> DensN() { return densn_; }

  //! routine to return density field at time step n-1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> DensNm() { return densnm_; }

  //! routine to return time deriv. of density at time step n for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> DensDtN() { return densdtn_; }

  //! routine to return time deriv. of density at time step n+1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> DensDtNp() { return densdtnp_; }

  //! routine to return discretization
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  //! return current time value
  double Time() const { return time_; }

  //! return current step number
  int Step() const { return step_; }

  //! return system matrix down-casted as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }

  //! return system matrix downcasted as block sparse matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }

  //! calculate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

protected:

  // don't want = operator and cctor
  //ScaTraTimIntImpl operator = (const ScaTraTimIntImpl& old);
  //! copy constructor
  ScaTraTimIntImpl (const ScaTraTimIntImpl& old);

  //! helper function to get algorithm title
  std::string MapTimIntEnumToString
  (
    const enum INPAR::SCATRA::TimeIntegrationScheme term  //!< the enum
  );

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  //! initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep() = 0;

  //! Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  //! reset the residual vector and add actual Neumann loads with time factor
  virtual void AddNeumannToResidual() = 0;

  //! add parameters specific for time-integration scheme
  virtual void AddSpecificTimeIntegrationParameters(ParameterList& params) = 0;

  //! Apply Dirichlet boundary conditions on provided state vector
  void ApplyDirichletBC
  (
    const double& time,  //!< evaluation time
    Teuchos::RCP<Epetra_Vector> phinp, //!< transported scalar(s) (may be = null)
    Teuchos::RCP<Epetra_Vector> phidt  //!< first time derivative (may be = null)
  );

  //! Apply Neumann boundary conditions
  void ApplyNeumannBC
  (
    const double& time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector>& phinp,   //!< transported scalar(s) state
    Teuchos::RCP<Epetra_Vector>& neumann_loads //!< Neumann loads
  );

  //! write state vectors (phinp and convective velocity) to BINIO
  void OutputState();

  //! write additional data required for restart
  virtual void OutputRestart() = 0;

  //! write flux vectors to BINIO
  void OutputFlux();

  //! routine to calculate mass/heat flux vectors
  Teuchos::RCP<Epetra_MultiVector> CalcFlux();

  /// the scalar transport discretization
  RCP<DRT::Discretization>   discret_;
  /// the solver
  RCP<LINALG::Solver>        solver_;
  /// the parameter list
  RCP<ParameterList>         params_;
  /// the discretization writer
  RCP<IO::DiscretizationWriter>  output_;
  /// processor id
  int                        myrank_;

  /// actual time
  double time_;

  /// actual step number
  int step_;

  /// flag defining type of scalar transport problem
  string prbtype_;

  /// maximum number of steps
  int stepmax_;

  /// maximum simulation time
  double maxtime_;

  /// used time integration scheme
  INPAR::SCATRA::TimeIntegrationScheme timealgo_;

  /// write results every upres_ steps
  int upres_;

  /// write restart data every uprestart_ steps
  int uprestart_;

  /// flag for writing the flux vector fields
  string writeflux_;

  /// time step size
  double dta_;

  /// time step size
  double dtp_;

  int cdvel_;

  /// maximum number of nonlinear iterations
  int itemax_;

  /// form of convective term
  string convform_;

  /// fine-scale subgrid-viscosity flag
  string fssgd_;

  // ELCH-specific paramter F/RT (default value for the temperature is 298K)
  double frt_;

  /// system matrix (either sparse matrix or block sparse matrix)
  RCP<LINALG::SparseOperator> sysmat_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  // maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// The residual vector
  RCP<Epetra_Vector>    residual_;

  /// phi at time n+1 and n
  RCP<Epetra_Vector>    phinp_;
  RCP<Epetra_Vector>    phin_;

  /// density at times n+1, n and n-1
  RCP<Epetra_Vector>    densnp_;
  RCP<Epetra_Vector>    densn_;
  RCP<Epetra_Vector>    densnm_;

  /// time derivative of density at times n+1 and n
  RCP<Epetra_Vector>  densdtnp_;
  RCP<Epetra_Vector>  densdtn_;

  /// density increment at time n+1
  RCP<Epetra_Vector>    densincnp_;

  /// histvector --- a linear combination of phinm, phin (BDF)
  ///                or phin, phidtn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  /// actual convective velocity
  RCP<Epetra_MultiVector>    convel_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    increment_;

  /// vm3 and avms solver
  RCP<FLD::VM3_Solver> vm3_solver_;
  RCP<FLD::AVMS_Solver> avms_solver_;

  /// only necessary for vm3 solver: fine-scale subgrid-diffusivity matrix
  RCP<LINALG::SparseMatrix> sysmat_sd_;

  /// only necessary for vm3 solver: subgrid-diffusivity-scaling vector
  RCP<Epetra_Vector>    subgrdiff_;

  //! Extractor used for convergence check of nonlinear iteration
  LINALG::MapExtractor      conpotsplitter_;

  /// Periodic Boundary Conditions. Guess what.
  RCP<PeriodicBoundaryConditions> pbc_;

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  // time measurement
  // ----------------
  double dtele_;
  double dtsolve_;

}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
#endif  // #ifdef CCADISCRET
