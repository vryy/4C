/*---------------------------------------------------------------------*/
/*!

\brief Optimization algorithm Globally Convergent Method of Moving Asyptotes (GCMMA)

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef OPTI_GCMMA_H_
#define OPTI_GCMMA_H_


#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_topopt.H"

class Epetra_MultiVector;
class Epetra_SerialDenseVector;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace IO

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
    explicit GCMMA(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> x, int numConstraints, Teuchos::RCP<Epetra_Vector> x_min,
        Teuchos::RCP<Epetra_Vector> x_max, Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// destructor
    virtual ~GCMMA(){};

    /// Perform one outer iteration step
    Teuchos::RCP<Epetra_Vector> Iterate(double& objective,
        Teuchos::RCP<Epetra_Vector> objectivegrad,
        Teuchos::RCP<Epetra_SerialDenseVector> constraints,
        Teuchos::RCP<Epetra_MultiVector> constraintsgrad);

    /// check if converged
    bool Converged(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
        Teuchos::RCP<Epetra_SerialDenseVector> constraints,
        Teuchos::RCP<Epetra_MultiVector> constraintsgrad);

    /// finish one iteration
    void FinishIteration(
        double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints, bool& doGradient);

    /// write output
    void Output();

    /// perform restart
    Teuchos::RCP<IO::DiscretizationReader> ReadRestart(int step);

    const Teuchos::RCP<const DRT::Discretization> Discretization() const { return discret_; }
    const Teuchos::RCP<IO::DiscretizationWriter> OutputWriter() const { return output_; }

    const Teuchos::RCP<const Epetra_Vector> X() const { return x_mma_; }
    Teuchos::RCP<Epetra_Vector> X() { return x_mma_; }
    double Obj() const { return obj_; }
    double Obj() { return obj_; }
    const Teuchos::RCP<const Epetra_Vector> ObjDeriv() const { return obj_deriv_; }
    Teuchos::RCP<Epetra_Vector> ObjDeriv() { return obj_deriv_; }
    const Teuchos::RCP<const Epetra_SerialDenseVector> Constr() const { return constr_; }
    Teuchos::RCP<Epetra_SerialDenseVector> Constr() { return constr_; }
    const Teuchos::RCP<const Epetra_MultiVector> ConstrDeriv() const { return constr_deriv_; }
    Teuchos::RCP<Epetra_MultiVector> ConstrDeriv() { return constr_deriv_; }

    int Iter() const { return total_iter_; }
    int OuterIter() const { return outer_iter_; }
    int InnerIter() const { return inner_iter_; }


   private:
    /// initialize one iteration
    void InitIter(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
        Teuchos::RCP<Epetra_SerialDenseVector> constraints,
        Teuchos::RCP<Epetra_MultiVector> constraintsgrad);

    /// update the moving asymptotes
    void Asymptotes();

    /// initialize rho in new outer iteration
    void InitRho();

    /// initialize rho in new inner iteration
    void UpdateRho(double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints);

    /// check convergence of inner iteration loop
    bool InnerConvergence(
        double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints, int& numNotFinished);

    /// initialize solution process of GCMMA subproblem
    void InitSubSolve();

    /// solution process of GCMMA subproblem
    void SubSolve();

    /// compute residual of GCMMA subproblem
    void ResApp(double& resnorm, double& resinf, const double& tol_sub);

    /// compute L2-norm of residuals
    double Res2Norm(Epetra_Vector* res1, Epetra_Vector* res2, Epetra_Vector* res3,
        Epetra_SerialDenseVector* res4, Epetra_SerialDenseVector* res5,
        Epetra_SerialDenseVector* res6, Epetra_SerialDenseVector* res7, double* res8, double* res9);

    /// compute LInf-norm of residuals
    double ResInfNorm(Epetra_Vector* res1, Epetra_Vector* res2, Epetra_Vector* res3,
        Epetra_SerialDenseVector* res4, Epetra_SerialDenseVector* res5,
        Epetra_SerialDenseVector* res6, Epetra_SerialDenseVector* res7, double* res8, double* res9);

    /// update variables after subsolve finished
    void Update();

    /// check KKT conditions for convergence
    bool KKTCond(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
        Teuchos::RCP<Epetra_SerialDenseVector> constraints,
        Teuchos::RCP<Epetra_MultiVector> constraintsgrad);

    void OutputToGmsh();


    //! @name problem data

    Teuchos::RCP<DRT::Discretization> discret_;  /// dof-map of optimization variable
    const Teuchos::ParameterList& params_;       /// optimizer parameter list
    //@}

    //! @name iteration counter

    int total_iter_;  /// counter for total iterations
    int outer_iter_;  /// counter for outer iterations
    int inner_iter_;  /// counter for inner iterations

    const int max_total_iter_;      /// maximal total iterations (outer and inner)
    const int max_inner_iter_;      /// maximal inner iterations
    const int max_outer_iter_;      /// maximal outer iterations (containing the gradient)
    const int max_sub_iter_;        /// maximal iterations within subproblem
    const int max_inner_sub_iter_;  /// maximal iterations of inner subproblem routine
    //@}

    //! @name size of optimization problem

    const int m_;      /// number of constraints
    const int n_loc_;  /// number of local optimization variables TODO auch bei periodischen RBs?
    const int
        n_;  /// number of optimization variables (=number of nodes) TODO auch bei periodischen RBs?
    //@}

    //! @name optimization variable data
    const INPAR::TOPOPT::DensityField dens_type_;  /// type of optimization field
    const INPAR::TOPOPT::SolverType solver_sub_;   /// solver type for solution of subproblem

    Teuchos::RCP<Epetra_Vector> x_;       /// current optimization variable
    Teuchos::RCP<Epetra_Vector> x_old_;   /// last optimization variable
    Teuchos::RCP<Epetra_Vector> x_old2_;  /// second to last optimization variable
    Teuchos::RCP<Epetra_Vector> x_mma_;   /// optimization variable of inner iteration

    Teuchos::RCP<Epetra_Vector> x_min_;  /// lower bound for objective variables
    Teuchos::RCP<Epetra_Vector> x_max_;  /// upper bound for objective variables

    Teuchos::RCP<Epetra_Vector> x_diff_;  /// bound difference for objective variables
    //@}

    //! @name moving asymptotes

    Teuchos::RCP<Epetra_Vector> asymp_min_;  /// lower moving asymptote
    Teuchos::RCP<Epetra_Vector> asymp_max_;  /// upper moving asymptote

    Teuchos::RCP<Epetra_Vector> alpha_;  /// lower bound for subproblem
    Teuchos::RCP<Epetra_Vector> beta_;   /// lower bound for subproblem
    //@}

    //! @name objective data

    double obj_;  /// objective value
    Teuchos::RCP<Epetra_Vector>
        obj_deriv_;    /// derivation of objective w.r.t. optimization variable
    double obj_appr_;  /// MMA approximation of objective
    //@}

    //! @name constraint data

    Teuchos::RCP<Epetra_SerialDenseVector> constr_;  /// constraint values
    Teuchos::RCP<Epetra_MultiVector>
        constr_deriv_;  /// derivation of constraints w.r.t. optimization variable
    Teuchos::RCP<Epetra_SerialDenseVector> constr_appr_;  /// MMA approximation of constraint
    //@}

    //! @name approximation function data

    Teuchos::RCP<Epetra_Vector> p0_;  /// first part of approximation function of objective
    Teuchos::RCP<Epetra_Vector> q0_;  /// second part of approximation function of objective
    double r0_;                       /// third part of approximation function of objective

    Teuchos::RCP<Epetra_MultiVector> P_;  /// first part of approximation function of constraints
    Teuchos::RCP<Epetra_MultiVector> Q_;  /// second part of approximation function of constraints
    Teuchos::RCP<Epetra_SerialDenseVector>
        b_;  /// third part of approximation function of constraints
    //@}

    //! @name approximation function coefficients

    double rho0_;  /// part of the objective function of the MMA subproblem
    Teuchos::RCP<Epetra_SerialDenseVector>
        rho_;               /// part of the objective function of the MMA subproblem, size m_
    const double rho0min_;  /// part of the objective function of the MMA subproblem
    Teuchos::RCP<Epetra_SerialDenseVector>
        rhomin_;             /// part of the objective function of the MMA subproblem
    const double rho_fac1_;  /// factor for updating rho
    const double rho_fac2_;  /// factor for updating rho
    //@}

    //! @name additional optimization variables

    Teuchos::RCP<Epetra_SerialDenseVector>
        y_mma_;     /// additional optimization variables for better convergence properties
    double z_mma_;  /// additional optimization variable for better convergence properties
    //@}

    //! @name Lagrange multiplier

    Teuchos::RCP<Epetra_Vector> xsi_;             /// Lagrange multiplier for xmin-x <= 0
    Teuchos::RCP<Epetra_Vector> eta_;             /// Lagrange multiplier for x-xmax <= 0
    Teuchos::RCP<Epetra_SerialDenseVector> lam_;  /// Lagrange multiplier for the m constraints
    Teuchos::RCP<Epetra_SerialDenseVector> mu_;   /// Lagrange multiplier for constraints y <= 0
    double zet_;                                  /// Lagrange multiplier for -z <= 0
    //@}

    //! @name additional model parameters

    double a0_;                                 /// model parameter
    Teuchos::RCP<Epetra_SerialDenseVector> a_;  /// model parameter with size = m_
    Teuchos::RCP<Epetra_SerialDenseVector> c_;  /// model parameter with size = m_
    Teuchos::RCP<Epetra_SerialDenseVector> d_;  /// model parameter with size = m_
    //@}

    //! @name convergence stuff

    double tol_sub_;                 /// tolerance for MMA subproblem
    const bool tol_sub_adaptiv_;     /// flag for adaptiv tolerance for MMA subproblem
    const double tol_sub_quot_fac_;  /// quotient factor for adaption of MMA subproblem tolerance
    const double tol_sub_min_;       /// minimal tolerance for MMA subproblem
    const double tol_kkt_;           /// tolerance for outer iterations (KKT-conditions are checked)
    const double fac_stepsize_;      /// factor for adjusting step size of optimization step
    const double asy_fac1_;          /// factor for computation of moving asymptotes
    const double asy_fac2_;          /// factor for computation of moving asymptotes
    const double fac_x_bd_;     /// factor for computation of boundaries for objective variables
    const double fac_sub_reg_;  /// regularisation factor in subproblem
    const double inc2_tol_;     /// factor for increasing u-x
    const double gamma_up_;     /// factor for increasing u-x
    const double gamma_down_;   /// factor for decreasing u-x
    //@}

    const double facmin_;                       /// minimal coefficient in a denominator
    Teuchos::RCP<Epetra_SerialDenseVector> s_;  /// helper variable

    //! @name output handling
    const int upres_;      /// write solution every upres steps to gmsh
    const int uprestart_;  /// write restart every uprestart steps

    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

  };  // end class GCMMA
}  // end namespace OPTI


#endif /* OPTI_GCMMA_H_ */
