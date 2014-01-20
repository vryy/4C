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

namespace ACOU
{
/*====================================================================*/
/*!
* Inverse Analysis for Photoacoustics
*/

class InvAnalysis
{
private:

  /// name given in command line
  std::string name_;

  /// scatra and acoustical discretizations to work with
  Teuchos::RCP<DRT::Discretization> scatra_discret_;
  Teuchos::RCP<DRT::DiscretizationHDG> acou_discret_;

  /// scatra and acoustical base algorithms to work with
  Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo_;
  Teuchos::RCP<ACOU::AcouImplicitTimeInt> acoualgo_;

  /// some things we need for the scatra and acou time integrators
  const Teuchos::ParameterList& scatraparams_;
  const Teuchos::ParameterList& scatrasolverparams_;
  Teuchos::RCP<Teuchos::ParameterList> acouparams_;
  Teuchos::RCP<LINALG::Solver> acousolver_;
  Teuchos::RCP<IO::DiscretizationWriter> acououtput_;
  Teuchos::RCP<IO::DiscretizationWriter> scatra_output_;

  /// time integration scheme
  INPAR::ACOU::DynamicType dyna_;

  /// processor id
  int myrank_;

  /// for iteration loop progress
  double error_;       /// actual error
  double tol_;         /// prescribed error tolerance
  double tol_grad_;    /// prescribed norm of gradient tolerance
  int iter_;           /// actual iteration number
  int max_iter_;       /// maximal iteration number for inverse analysis
  int max_ls_iter_;    /// maximal iteration number for line search
  int output_count_;   /// auxiliary variable for nicely named output files
  bool fdcheck_;       /// flog indicating performance of finite difference check

  /// regularization coefficient for absorption coefficient
  double alpha_;

  /// flag whether we want to do check for the gradient as convergence criteria
  bool calcacougrad_;

  /// container for right hand side values of acoustical adjoint problem
  Teuchos::RCP<Epetra_MultiVector> acou_rhs_;
  Teuchos::RCP<Epetra_MultiVector> acou_rhsm_;
  Teuchos::RCP<Epetra_Map> abcnodes_map_;

  /// map extractor for monitored / abc nodes
  Teuchos::RCP<LINALG::MapExtractor> abcnodes_mapex_;

  /// monitor file related things
  unsigned int         nnodes_;        /// number of nodes observed
  int                  nsteps_;        /// number of time steps observed
  std::vector<int>     nodes_;         /// list of node gids observed
  std::vector<double>  timesteps_;     /// time values where pressure values are measured
  Epetra_SerialDenseVector mcurve_;    /// measured pressure values of the experiments (target value)

  /// solution of the acoustical adjoint problem at initial step (so to say last step of the backward time integration)
  Teuchos::RCP<Epetra_Vector> adjoint_phi_0_;

  /// solution of the optical forward problem
  Teuchos::RCP<Epetra_Vector> phi_;
  /// solution of the optical adjoint problem
  Teuchos::RCP<Epetra_Vector> adjoint_w_;

  /// node based vector holding material properties
  Teuchos::RCP<Epetra_Vector> node_c_;     /// speed of sound
  Teuchos::RCP<Epetra_Vector> node_rho_;   /// density
  Teuchos::RCP<Epetra_Vector> node_mu_a_;  /// absorption coefficient


  Epetra_SerialDenseVector p_;        /// parameter vector
  Epetra_SerialDenseVector pm_;       /// parameter vector from last iteration
  int np_;                            /// number of parameters
  int nm_;                            /// number of materials
  double t_;                          /// maximal time of acoustical integration
  double dtacou_;                     /// time step size for acoustical simulations


  Epetra_SerialDenseVector G_;        /// the gradient
  Epetra_SerialDenseVector Gm_;       /// gradient last iteration

  /// the approximation of the inverse Hessian (approximated via BFGS)
  Epetra_SerialDenseMatrix H_;
  Epetra_SerialDenseMatrix Hm_;       /// approximated Hessian from last iteration

  /// the objective function value
  double J_;

  /// norm of the difference of old and new parameters (which can also be used as stopping criteria)
  double normdiffp_;

public:

  /// standard Constructor
  InvAnalysis(
    Teuchos::RCP<DRT::Discretization> scatradis, ///< the light discretization
    Teuchos::RCP<DRT::DiscretizationHDG> acoudis, ///< the sound discretization
    const Teuchos::ParameterList& scapara,
    const Teuchos::ParameterList& scasolvpara,
    Teuchos::RCP<Teuchos::ParameterList> acoupara,
    Teuchos::RCP<LINALG::Solver> acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  /// Read in routine for the material parameters
  void ReadInParameters();

  /// Set material parameters
  void SetParameters(Epetra_SerialDenseVector p_cur);

  /// Actual inverse analysis loop
  void Integrate();

  /// create test field
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  /// for result test frame work
  Teuchos::RCP<DRT::Discretization> ScatraDiscretization() { return scatra_discret_;};

  /// output of element based absorption coefficient vector for resulttest
  const Teuchos::RCP<Epetra_Vector> ElementMatVec();

  /// The underlying subproblems
  //@{
  /// Solve the standard problem
  void SolveStandardProblem();

  /// Calculate the value of the objective function (after solving the forward problem)
  void CalculateObjectiveFunctionValue();

  /// Solve the acoustical adjoint problem
  void SolveAdjointAcouProblem();

  /// Solve the optical adjoint problem
  void SolveAdjointOptiProblem();

  /// Build the gradient
  void CalculateGradient();

  /// Finite difference check on gradient
  void FD_GradientCheck();

  /// Update the BFGS approximation of the Hessian
  void UpdateHessian();

  /// Update the sought parameters
  void UpdateParameters();

  /// Line search routine to solve the internal BFGS problem
  Epetra_SerialDenseVector LineSearch(Epetra_SerialDenseVector d);

  /// Tell all the elements their new material parameters
  void UpdateMaterial();

  /// Calculate the error and the value of the objective function, update iteration count, update Gradient
  void CalculateStatsAndService();

  /// Output some useful user information, like time consume, solution advance, ...
  void OutputStats();

  /// Check for convergence, returns false when converged
  bool ConvergenceCheck();
  //@}

}; // class InvAnalysis

} // namespace ACOU



#endif /* ACOU_INV_ANALYSIS_H */
