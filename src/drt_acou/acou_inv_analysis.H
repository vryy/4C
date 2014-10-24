/*----------------------------------------------------------------------*/
/*!
 * \file acou_inv_analysis.H

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ACOU_INV_ANALYSIS_H
#define ACOU_INV_ANALYSIS_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_inpar/inpar_acou.H"

// forward declarations

namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class ResultTest;
}
namespace ACOU
{
  class AcouImplicitTimeInt;
  class PatMatParManager;
}
namespace SCATRA
{
  class TimIntStationary;
}
namespace LINALG
{
  class Solver;
  class MapExtractor;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace INVANA
{
  class MatParManager;
}


namespace DRT
{
  namespace UTILS
  {
    template <typename> class TimIntMStep;
  }
}

namespace ACOU
{
/*====================================================================*/
/*!
* Inverse Analysis for Photoacoustics
*/

class InvAnalysis
{

public:

  //! standard Constructor
  InvAnalysis(
    Teuchos::RCP<DRT::Discretization> scatradis, //!< the light discretization
    Teuchos::RCP<DRT::DiscretizationHDG> acoudis, //!< the sound discretization
    const Teuchos::ParameterList& scapara,
    const Teuchos::ParameterList& scasolvpara,
    Teuchos::RCP<Teuchos::ParameterList> acoupara,
    Teuchos::RCP<LINALG::Solver> acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! Actual inverse analysis loop
  void Integrate();

  //! create test field
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  //! for result test frame work
  Teuchos::RCP<DRT::Discretization> ScatraDiscretization() { return scatra_discret_; }

  //! output of element based absorption coefficient vector for resulttest
  const Teuchos::RCP<Epetra_MultiVector> ElementMatVec();

private:

  //! name given in command line
  std::string name_;

  //! scatra and acoustical discretizations to work with
  Teuchos::RCP<DRT::Discretization> scatra_discret_;
  Teuchos::RCP<DRT::DiscretizationHDG> acou_discret_;

  //! scatra and acoustical base algorithms to work with
  Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo_;
  Teuchos::RCP<ACOU::AcouImplicitTimeInt> acoualgo_;

  //! some things we need for the scatra and acou time integrators
  const Teuchos::ParameterList& scatraparams_;
  const Teuchos::ParameterList& scatrasolverparams_;
  Teuchos::RCP<Teuchos::ParameterList> acouparams_;
  Teuchos::RCP<LINALG::Solver> acousolver_;
  Teuchos::RCP<IO::DiscretizationWriter> acououtput_;
  Teuchos::RCP<IO::DiscretizationWriter> scatra_output_;

  //! time integration scheme
  INPAR::ACOU::DynamicType dyna_;
  //! physical type
  INPAR::ACOU::PhysicalType phys_;
  //! optimization tyoe
  INPAR::ACOU::OptimizationType opti_;

  //! processor id
  int myrank_;

  //! for iteration loop progress
  double error_;       //! actual error
  double tol_;         //! prescribed error tolerance
  double tol_grad_;    //! prescribed norm of gradient tolerance
  int iter_;           //! actual iteration number
  int max_iter_;       //! maximal iteration number for inverse analysis
  int max_ls_iter_;    //! maximal iteration number for line search
  int output_count_;   //! auxiliary variable for nicely named output files
  bool fdcheck_;       //! flag indicating performance of finite difference check
  double ls_c_;        //! parameter for calculation of sufficient decrease condition
  double ls_rho_;      //! step length reduction given by this parameter
  double ls_gd_scal_;  //! factor for line search step length scaling

  //! regularization coefficient for absorption coefficient
  double alpha_;
  double beta_;

  //! flag whether we want to do check for the gradient as convergence criteria
  bool calcacougrad_;

  //! gradient scaling
  bool scalegradele_;

  //! container for right hand side values of acoustical adjoint problem
  Teuchos::RCP<Epetra_MultiVector> acou_rhs_;
  Teuchos::RCP<Epetra_MultiVector> acou_rhsm_;
  Teuchos::RCP<Epetra_Map> abcnodes_map_;

  //! map extractor for monitored / abc nodes
  Teuchos::RCP<LINALG::MapExtractor> abcnodes_mapex_;

  //! list of observed node gids
  std::vector<int>     nodes_;

  //! solution of the acoustical adjoint problem at initial step (so to say last step of the backward time integration)
  Teuchos::RCP<Epetra_Vector> adjoint_phi_0_;
  //! solution of the optical forward problem
  Teuchos::RCP<Epetra_Vector> phi_;
  //! solution of the optical adjoint problem
  Teuchos::RCP<Epetra_Vector> adjoint_w_;

  //! material handling done by the MatParManager from Sebastian Kehl
  Teuchos::RCP<INVANA::MatParManager> matman_;

  //! gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_;

  //! previous gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_o_;

  //! current direction of steepest descent maybe properly scaled by hessian approximation
  Teuchos::RCP<Epetra_MultiVector> direction_;

  //! current step based on direction direction_ scaled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> step_;

  //! current direction of steepest descent maybe properly scaled by hessian approximation
  Teuchos::RCP<Epetra_MultiVector> d_;

  //! norm of the gradient at current iteration
  double normgrad_;

  //! norm of the gradient in first iteration
  double normgrad0_;

  //! maximum size of storage (number of multivectors)
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! storage of solution delta
  Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > sstore_;

  //! storage of gradient delta
  Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > ystore_;

  //! node based vector holding absorption coefficient
  Teuchos::RCP<Epetra_Vector> node_mu_a_;

  //! time step size for acoustical simulations
  double dtacou_;

  //! the objective function value
  double J_;

  //! norm of the difference of old and new parameters (which can also be used as stopping criteria)
  double normdiffp_;

  //! time measurement
  double tstart_;

  //! Read in measured values from the monitor file
  void ReadMonitor(std::string monitorfilename);

  //! Compute node based vectors
  void ComputeNodeBasedVectors();

  //! The underlying subproblems
  //@{
  //! Solve the standard problem
  void SolveStandardProblem();

  //! Calculate the value of the objective function (after solving the forward problem)
  void CalculateObjectiveFunctionValue();

  //! Solve the acoustical adjoint problem
  void SolveAdjointAcouProblem();

  //! Solve the optical adjoint problem
  void SolveAdjointOptiProblem();

  //! Calculate RhsVec for adjoint Opti Problem
  Teuchos::RCP<Epetra_Vector> CalculateNonconformRhsvec(Teuchos::RCP<Epetra_Vector> acounodevec);

  //! Build the gradient
  void CalculateGradient();

  //! Finite difference check on gradient
  void FD_GradientCheck();

  //! store new vector quantities
  void StoreVectors();

  //! compute new direction based on lbfgs formula
  void ComputeDirection();

  //! Update the sought parameters, return argument is true if successful
  bool UpdateParameters();

  //! Line search routine to solve the internal BFGS problem, return argument is true if successful
  bool LineSearch();

  //! Calculate the error and the value of the objective function, update iteration count, update Gradient
  void CalculateStatsAndService();

  //! Output some useful user information, like time consume, solution advance, ...
  void OutputStats();

  //! Check for convergence, returns false when converged
  bool ConvergenceCheck();
  //@}

}; // class InvAnalysis

} // namespace ACOU



#endif /* ACOU_INV_ANALYSIS_H */
