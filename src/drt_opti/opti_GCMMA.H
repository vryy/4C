/*!------------------------------------------------------------------------------------------------*
\file opti_GCMMA.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef OPTI_GCMMA_H_
#define OPTI_GCMMA_H_


#include<Teuchos_RCP.hpp>

class Epetra_MultiVector;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace Teuchos
{
  class ParameterList;
}


namespace OPTI
{
  /*
  \brief Optimization approach based on the globally convergent method of moving asymptotes

  This algorithm is well suited for optimization problems with high dimensional
  optimization variables (for example discretized space-dependent functions)

  The number of constraints is assumed to be moderate, but a high number of
  constraints should also be handled efficiently.

  The method of moving asymptotes belongs to the class of optimization approaches
  based on conservative convex separable approximations. A series of convex
  subproblems is solved (which can be done efficiently) so that in every iteration
  the objective becomes lower.

  The globally convergent version GCMMA of the MMA performs additional outer
  iterations which results in global convergence of the method.

  See Svanberg 2002, "A class of globally convergent optimization methods based
  on conservative convex separable approximations" for details.

  \author winklmaier
  \date 12/11
  */
  class GCMMA
  {
  public:
    /// constructor
    explicit GCMMA(
        Teuchos::RCP<const DRT::Discretization> discret,
        const Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> x,
        int numConstraints,
        Teuchos::RCP<Epetra_Vector> x_min,
        Teuchos::RCP<Epetra_Vector> x_max
    );

    /// desctructor
    virtual ~GCMMA(){};

    /// Perform one outer iteration step
    void Iterate(
        double& objective,
        Teuchos::RCP<Epetra_Vector> objectivegrad,
        double* constraints,
        Teuchos::RCP<Epetra_MultiVector> constraintsgrad,
        bool& innerIterConverged);

    bool Converged();

  private:

    void PrepareOuterIter();

    void PrepareInnerIter();

    void Asymptotes();

    void UpdateRho();

    bool InnerConvergence();

    void InitSubSolve();

    void SubSolve();

    void Update(
        bool& innerIterConverged
    );

    //! @name problem data

    Teuchos::RCP<const DRT::Discretization> discret_; /// dof-map of optimization variable
    const Teuchos::ParameterList& params_; /// optimization parameter list
    //@}

    //! @name iteration counter

    int outer_iter_; /// counter for outer iterations
    int inner_iter_; /// counter for inner iterations

    const int max_outer_iter_; /// maximal outer iterations
    const int max_inner_iter_; /// maximal inner iterations
    //@}

    //! @name size of optimization problem

    const int m_; /// number of constraints
    const int n_; /// number of optimization variables (=number of nodes) TODO auch bei periodischen RBs
    //@}

    //! @name optimization variable data

    Teuchos::RCP<Epetra_Vector> x_; /// current optimization variable
    Teuchos::RCP<Epetra_Vector> x_old_; /// last optimization variable
    Teuchos::RCP<Epetra_Vector> x_old2_; /// second to last optimization variable
    Teuchos::RCP<Epetra_Vector> x_mma_; /// optimization variable of inner iteration

    Teuchos::RCP<Epetra_Vector> x_min_; /// lower bound for objective variables
    Teuchos::RCP<Epetra_Vector> x_max_; /// upper bound for objective variables

    Teuchos::RCP<Epetra_Vector> x_diff_; /// bound difference for objective variables
    double x_diff_min_; /// minimal difference of x_min and x_max for divisions
    //@}

    //! @name moving asymptotes

    Teuchos::RCP<Epetra_Vector> asymp_min_; /// lower moving asymptote
    Teuchos::RCP<Epetra_Vector> asymp_max_; /// upper moving asymptote

    Teuchos::RCP<Epetra_Vector> alpha_; /// lower bound for subproblem
    Teuchos::RCP<Epetra_Vector> beta_; /// lower bound for subproblem
    //@}

    //! @name objective data

    double obj_; /// objective value
    double obj_mma_; /// objective value of inner iteration

    Teuchos::RCP<Epetra_Vector> obj_deriv_; /// derivation of objective w.r.t. optimization variable
    double obj_appr_; /// MMA approximation of objective
    //@}

    //! @name constraint data

    double* constr_; /// constraint values
    double* constr_mma_; /// constraint value of inner iteration

    Teuchos::RCP<Epetra_MultiVector> constr_deriv_; /// derivation of constraints w.r.t. optimization variable
    double* constr_appr_; /// MMA approximation of constraint
    //@}

    //! @name approximation function data

    Teuchos::RCP<Epetra_Vector> p0_; /// first part of approximation function of objective
    Teuchos::RCP<Epetra_Vector> q0_; /// second part of approximation function of objective
    double r0_; /// third part of approximation function of objective

    Teuchos::RCP<Epetra_MultiVector> P_; /// first part of approximation function of constraints
    Teuchos::RCP<Epetra_MultiVector> Q_; /// second part of approximation function of constraints
    double* r_; /// third part of approximation function of constraints
    //@}

    //! @name approximation function coefficients

    double rho0_; /// part of the objective function of the MMA subproblem
    double* rho_; /// part of the objective function of the MMA subproblem, size m_
    double rho0min_; /// part of the objective function of the MMA subproblem
    double* rhomin_; /// part of the objective function of the MMA subproblem
    //@}

    //! @name additional optimization variables

    double* y_mma_; /// additional optimization variables for better convergence properties
    double z_mma_; /// additional optimization variable for better convergence properties
    //@}

    //! @name Lagrange multiplier

    Teuchos::RCP<Epetra_Vector> xsi_; /// Lagrange multiplier for xmin-x <= 0
    Teuchos::RCP<Epetra_Vector> eta_; /// Lagrange multiplier for x-xmax <= 0
    double* lam_; /// Lagrange multiplier for the m constraints
    double* mu_; /// Lagrange multiplier for constraints y <= 0
    double zet_; /// Lagrange multiplier for -z <= 0
    //@}

    //! @name additional model parameters

    double a0_; /// model parameter
    double* a_; /// model parameter with size = m_
    double* c_; /// model parameter with size = m_
    double* d_; /// model parameter with size = m_
    //@}

    //! @name convergence stuff

    double tol_sub_min_; /// tolerance for MMA subproblem

    double tol_; /// tolerance for outer iterations (KKT-conditions are checked)
    double res_norm_; /// norm of the residual (KKT-conditions norm)
    double inc_norm_; /// norm of the increment
    //@}


    double facmin_; /// minimal coefficient in a denominator
    double* s_; /// helper variable



  }; // end class GCMMA
} // end namespace OPTI


#endif /* OPTI_GCMMA_H_ */
