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
  class AcouTimeInt;
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
  class RegularizationBase;
}

namespace TIMINT
{
  template <typename> class TimIntMStep;
}

namespace ACOU
{

class PATLineSearch;

/*====================================================================*/
/*!
* Inverse Analysis for Photoacoustics
*/
class InvAnalysis
{
  friend class PATLineSearch;

public:

  //! standard Constructor
  InvAnalysis(
    Teuchos::RCP<DRT::Discretization> scatradis, //!< the light discretization
    Teuchos::RCP<DRT::DiscretizationHDG> acoudis, //!< the sound discretization
    Teuchos::RCP<Teuchos::ParameterList> scatrapara,
    Teuchos::RCP<LINALG::Solver> scatrasolv,
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
  Teuchos::RCP<DRT::Discretization>    scatra_discret_;
  Teuchos::RCP<DRT::DiscretizationHDG> acou_discret_;

  //! scatra and acoustical base algorithms to work with
  Teuchos::RCP<SCATRA::TimIntStationary>  scatraalgo_;
  Teuchos::RCP<ACOU::AcouTimeInt> acoualgo_;

  //! some things we need for the scatra and acou time integrators
  Teuchos::RCP<Teuchos::ParameterList>   scatraparams_;
  Teuchos::RCP<Teuchos::ParameterList>   scatraextraparams_;
  Teuchos::RCP<LINALG::Solver>           scatrasolver_;
  Teuchos::RCP<Teuchos::ParameterList>   acouparams_;
  Teuchos::RCP<LINALG::Solver>           acousolver_;
  Teuchos::RCP<IO::DiscretizationWriter> acououtput_;
  Teuchos::RCP<IO::DiscretizationWriter> scatraoutput_;

  //! time integration scheme
  INPAR::ACOU::DynamicType      dyna_;
  //! physical type
  INPAR::ACOU::PhysicalType     phys_;
  //! optimization type
  INPAR::ACOU::OptimizationType opti_;

  //! processor id
  int myrank_;

  //! for iteration loop progress
  double error_;          //! actual error
  double error_start_;    //! initial error
  double tol_;            //! prescribed error tolerance
  double tol_grad_;       //! prescribed norm of gradient tolerance
  int    iter_;           //! actual iteration number
  int    max_iter_;       //! maximal iteration number for inverse analysis
  int    max_ls_iter_;    //! maximal iteration number for line search
  int    output_count_;   //! auxiliary variable for nicely named output files
  int    last_acou_out_;  //! last forward acoustic output count
  bool   fdcheck_;        //! flag indicating performance of finite difference check
  double ls_c_;           //! parameter for calculation of sufficient decrease condition
  double ls_rho_;         //! step length reduction given by this parameter
  double ls_gd_scal_;     //! factor for line search step length scaling
  bool   backprojection_; //! flag if material should be estimated with backprojection
  bool   meshconform_;    //! is acoustical mesh conform with optical

  //! flag whether we want to do check for the gradient as convergence criteria
  bool calcacougrad_;

  //! gradient scaling
  bool scalegradele_;

  //! flag indicating if acoustic material is also optimized
  bool acouopt_;

  //! bool indicating if segmentation is performed
  bool segmentation_;
  double penaltyparam_;

  //! contribution from regularization
  double adaptreguval_;
  double reguweight_;

  //! container for right hand side values of acoustical adjoint problem
  Teuchos::RCP<Epetra_MultiVector> acou_rhs_;
  Teuchos::RCP<Epetra_MultiVector> acou_rhsm_;
  Teuchos::RCP<Epetra_Map>         abcnodes_map_;

  //! map extractor for monitored / abc nodes
  Teuchos::RCP<LINALG::MapExtractor> abcnodes_mapex_;

  //! list of observed node gids
  std::vector<int>                  nodes_;
  Epetra_SerialDenseVector          nod_curve_interpol_;
  std::vector<double>               timesteps_;

  //! solution of the acoustical adjoint problem at initial step (so to say last step of the backward time integration)
  Teuchos::RCP<Epetra_Vector> adjoint_psi_;
  //! solution of the optical forward problem
  Teuchos::RCP<Epetra_Vector> phi_;
  //! solution of the optical adjoint problem
  Teuchos::RCP<Epetra_Vector> adjoint_phi_;

  //! material handling done by the MatParManager from Sebastian Kehl
  Teuchos::RCP<INVANA::MatParManager> scatra_matman_;
  //! material handler in case of acoustical variable materials
  Teuchos::RCP<INVANA::MatParManager> acou_matman_;

  //!regularization done by the RegularizationBase from Sebastian Kehl
  Teuchos::RCP<INVANA::RegularizationBase> scatra_regman_;
  Teuchos::RCP<INVANA::RegularizationBase> acou_regman_;

  //! gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> scatra_objgrad_;
  Teuchos::RCP<Epetra_MultiVector> acou_objgrad_;

  //! previous gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> scatra_objgrad_o_;
  Teuchos::RCP<Epetra_MultiVector> acou_objgrad_o_;

  //! current step based on direction direction_ scaled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> scatra_step_;
  Teuchos::RCP<Epetra_MultiVector> acou_step_;

  //! current direction of steepest descent maybe properly scaled by hessian approximation
  Teuchos::RCP<Epetra_MultiVector> scatra_d_;
  Teuchos::RCP<Epetra_MultiVector> acou_d_;

  //! vector holding the f values
  Teuchos::RCP<Epetra_MultiVector> k_;

  //! current direction of steepest descent for f
  Teuchos::RCP<Epetra_MultiVector> k_d_;

  //! current gradient for f
  Teuchos::RCP<Epetra_MultiVector> k_objgrad_;

  //! line search class
  Teuchos::RCP<PATLineSearch> linesearch_;

  //! norm of the gradient at current iteration
  double scatra_normgrad_;
  double acou_normgrad_;

  //! norm of the gradient in first iteration
  double scatra_normgrad0_;
  double acou_normgrad0_;

  //! maximum size of storage (number of multivectors)
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! storage of solution delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > scatra_sstore_;
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > acou_sstore_;

  //! storage of gradient delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > scatra_ystore_;
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > acou_ystore_;

  //! node based vector holding absorption coefficient
  Teuchos::RCP<Epetra_Vector> node_mu_a_;

  //! acoustical time step
  double dtacou_;

  //! the objective function value
  double J_;

  //! the initial objective function value
  double J_start_;

  //! norm of the difference of old and new parameters (which can also be used as stopping criteria)
  double normdiffp_;

  //! time measurement
  double tstart_;

  //! Read in measured values from the monitor file
  void ReadMonitor(std::string monitorfilename);

  //! Calculate the distance between mics and nods in x,y,z
  double Delta(double coord_M_x,double coord_M_y, double coord_M_z, double coord_N_x,double coord_N_y, double coord_N_z);

  //! Write the interpolated data
  void Interpol(const double nod_coords[3],std::vector<std::vector<double> > mic_coords, unsigned int mic_1, unsigned int mic_2, int nmic, unsigned int nod, int timesteps, Epetra_SerialDenseVector& curve, Epetra_SerialDenseVector& inter_curve);

  //! Calculate epsilon individually
  double GetEpsilon(int nnodes);

  //! Estimate the mu values of each element via back-projection
  void EstimateMua(int count);

  //! Compute node based vectors
  void ComputeNodeBasedVectors();

  //! The underlying subproblems
  //@{
  //! optical forward problem
  void SolveStandardProblemScatra();

  //! output for parameters
  void OutputReactionAndDiffusion();

  //! acoustical forward problem
  void SolveStandardProblemAcou();

  //! Calculate the value of the objective function (after solving the forward problem)
  void CalculateObjectiveFunctionValue(bool first);

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

  //! Calculate the error and the value of the objective function, update iteration count, update Gradient
  void CalculateStatsAndService();

  //! Output some useful user information, like time consume, solution advance, ...
  void OutputStats();

  //! Check for convergence, returns false when converged
  bool ConvergenceCheck();
  //@}

  struct Segmentation
  {
    // number of materials
    unsigned int nummats;

    // table holding the material values
    std::vector<std::vector<double> > materialtable;

    // penalty parameter
    double penalty;

    // pointer to the scatra and acoustic discretization
    Teuchos::RCP<DRT::Discretization> scadisc;
    Teuchos::RCP<DRT::DiscretizationHDG> acoudisc;

    // constructor
    Segmentation(std::string materialfilename,
                 double penaltypara,
                 Teuchos::RCP<DRT::Discretization> scatradiscretization,
                 Teuchos::RCP<DRT::DiscretizationHDG> acoudiscretization);

    // initialization of the parameters
    void InitParameters(Teuchos::RCP<Epetra_MultiVector> concentrations,
                        Teuchos::RCP<INVANA::MatParManager> scatramatman,
                        Teuchos::RCP<INVANA::MatParManager> acoumatman,
                        Teuchos::RCP<Epetra_MultiVector> scatraparams,
                        Teuchos::RCP<Epetra_MultiVector> acouparams);

    // get as input the current parameters and gradients
    // use all of them to calculate the new step direction
    // considering f
    void CalculateGradient(
        Teuchos::RCP<Epetra_MultiVector> sca_params,
        Teuchos::RCP<Epetra_MultiVector> acou_params,
        Teuchos::RCP<Epetra_MultiVector> sca_grad,
        Teuchos::RCP<Epetra_MultiVector> acou_grad,
        Teuchos::RCP<Epetra_MultiVector> f_params,
        Teuchos::RCP<Epetra_MultiVector> f_grad);

    // update of the parameters
    void UpdateParameters(Teuchos::RCP<Epetra_MultiVector> k_p,
                          Teuchos::RCP<INVANA::MatParManager> scatramatman,
                          Teuchos::RCP<INVANA::MatParManager> acoumatman,
                          Teuchos::RCP<Epetra_MultiVector> sca_p,
                          Teuchos::RCP<Epetra_MultiVector> acou_p);
  };

  Teuchos::RCP<Segmentation> segmentationsolver_;

}; // class InvAnalysis

/// Line search class specific to PAT
/*!
    This is a class providing line search functionality. We do not use the standard
    class from the nonlinear solver, since this one has to be suited to PAT needs, e.g.
    to choose between simultaneous parameter update or acoustic / optic one after another.

    This implementation is based on the description in the book "Numerical Optimization",
    J. Nocedal, J. Wrigth, Chapter 3.5.

    \author schoeder
*/
class PATLineSearch
{

public:

  //! standard constructor
  PATLineSearch(Teuchos::RCP<InvAnalysis> problem, bool successive);

  //! initialization
  void Init();

  //! actual line search routine
  bool Run();

private:

  //! pointer to the inverse problem
  Teuchos::RCP<InvAnalysis> invproblem_;

  //! maximum number of iterations
  int itermax_;

  //! maximal step length
  double alpha_max_;

  //! bool indicating if also acoustical parameters are optimized
  bool acouopt_;

  //! bool indicating if successive or simultaneous update is performed
  bool successive_;

  //! bool indicating if segmentation is performed
  bool segmentation_;

  //! setting at the beginning of an iteration
  double J_0_;        //! initial objective function value
  double normsca_0_;  //! initial norm of the scatra search (grad*dir)
  double normacou_0_; //! initial norm of the acou search (grad*dir)
  double normk_0_;    //! initial norm of f direction -(dir*dir)=grad*dir

  //! variables during an iteration
  double J_im1_;      //! previous objective function value
  double alpha_x_;    //! the final step length
  double alpha_i_;    //! the current step length
  double alpha_im1_;  //! the previous step length
  double normsca_;    //! norm of scatra search
  double normacou_;   //! norm of acou search
  double normk_;      //! norm of f search

  //! parameters for the line search
  double c1_;
  double c2_;

  //! zoom function
  double Zoom(double alpha_lo, double alpha_hi, double J_0, double J_alpha_lo, double normphi_0, bool turn, bool sca, bool acou);

  //! different run routines: successive or simultaneous
  bool RunSuccessive();
  bool RunSimultaneous();
  bool RunSegmentation();

  //! initial f values
  Teuchos::RCP<Epetra_MultiVector> k_o_;

}; // class PATLineSearch

} // namespace ACOU



#endif /* ACOU_INV_ANALYSIS_H */
