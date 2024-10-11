/*----------------------------------------------------------------------*/
/*! \file
\brief Class containing Uzawa algorithm to solve linear system.
\level 2
*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_SOLVER_HPP
#define FOUR_C_CONSTRAINT_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class MapExtractor;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace CONSTRAINTS
{
  /*!
  \brief Class containing uzawa algorithm to solve linear system.
  */
  class ConstraintSolver
  {
   public:
    /*!
    \brief Constructor
    */
    ConstraintSolver(Teuchos::RCP<Core::FE::Discretization> discr,  ///< discretization
        Core::LinAlg::Solver& solver,  ///< Solver to solve linear subproblem in iteration
        Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,  ///< Map extractor for Dirichlet DOFs
        Teuchos::ParameterList param  ///< parameterlist containing solver parameters
    );


    /*!
      \brief Set it up
     */
    void setup(Core::FE::Discretization& discr, Core::LinAlg::Solver& solver,
        Core::LinAlg::MapExtractor& dbcmaps, Teuchos::ParameterList params);

    /*!
      \brief Solve constraint linear system
    */
    void solve(Core::LinAlg::SparseMatrix& stiff,  ///< stiffness matrix
        Core::LinAlg::SparseMatrix& constr,        ///< constraint matrix with Dirichlet zeros
        Core::LinAlg::SparseMatrix&
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Core::LinAlg::Vector<double>> dispinc,  ///< displacement increment to compute
        Core::LinAlg::Vector<double>& lagrinc,      ///< lagrange multiplier increment to compute
        Core::LinAlg::Vector<double>& rhsstandard,  ///< standard right hand side
        Core::LinAlg::Vector<double>& rhsconstr     ///< constraint errors
    );

    /*!
       \brief Return the current value of the uzawa parameter
    */
    double get_uzawa_parameter() { return iterationparam_; }

    /*!
      \brief Set uzawa parameter (used for dynamic adaptation and restart)
    */
    void set_uzawa_parameter(double restartval  ///< value to replace Uzawa Parameter with
    )
    {
      iterationparam_ = restartval;
      counter_ = 1;
      return;
    }

    void set_stc_prop(
        Inpar::Solid::StcScale stcalgo, Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat)
    {
      stcalgo_ = stcalgo;
      stcmat_ = stcmat;
    };

   private:
    // do not want = operator, cctor
    ConstraintSolver operator=(const ConstraintSolver& old);
    ConstraintSolver(const ConstraintSolver& old);

    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;
    */
    void solve_uzawa(Core::LinAlg::SparseMatrix& stiff,  ///< stiffness matrix
        Core::LinAlg::SparseMatrix& constr,              ///< constraint matrix with Dirichlet zeros
        Core::LinAlg::SparseMatrix&
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Core::LinAlg::Vector<double>> dispinc,  ///< displacement increment to compute
        Core::LinAlg::Vector<double>& lagrinc,      ///< lagrange multiplier increment to compute
        Core::LinAlg::Vector<double>& rhsstandard,  ///< standard right hand side
        Core::LinAlg::Vector<double>& rhsconstr     ///< constraint errors
    );

    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;
    */
    void solve_simple(Core::LinAlg::SparseMatrix& stiff,  ///< stiffness matrix
        Core::LinAlg::SparseMatrix& constr,  ///< constraint matrix with Dirichlet zeros
        Core::LinAlg::SparseMatrix&
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Core::LinAlg::Vector<double>& dispinc,      ///< displacement increment to compute
        Core::LinAlg::Vector<double>& lagrinc,      ///< lagrange multiplier increment to compute
        Core::LinAlg::Vector<double>& rhsstandard,  ///< standard right hand side
        Core::LinAlg::Vector<double>& rhsconstr     ///< constraint errors
    );

    /*!
      \brief Solve linear system directly by assembling everything into one big matrix
    */
    void solve_direct(Core::LinAlg::SparseMatrix& stiff,  ///< stiffness matrix
        Core::LinAlg::SparseMatrix& constr,  ///< constraint matrix with Dirichlet zeros
        Core::LinAlg::SparseMatrix&
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Core::LinAlg::Vector<double>& dispinc,      ///< displacement increment to compute
        Core::LinAlg::Vector<double>& lagrinc,      ///< lagrange multiplier increment to compute
        Core::LinAlg::Vector<double>& rhsstandard,  ///< standard right hand side
        Core::LinAlg::Vector<double>& rhsconstr     ///< constraint errors
    );

    Teuchos::RCP<Core::FE::Discretization> actdisc_;  ///< standard discretization
    int max_iter_;                                    ///< number of maximal iterations
    double iterationparam_;                           ///< parameter for Uzawa algorithm
    double minparam_;      ///< minimal possible parameter for Uzawa algorithm
    double iterationtol_;  ///< tolerance
    double tolres_;        ///< tolerance for residual
    double tolconstr_;     ///< tolerance for constraint
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        dirichtoggle_;                                  ///< \b only for compatability: dirichlet
                                                        ///< toggle -- monitor its target change!
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        firstdispinc_;  ///< history variable holding displacement increment for first NRI
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        firstlagrinc_;    ///< history variable holding multiplier increment for first NRI
    bool isadapttol_;     ///< adaptive tolerance for solver?
    bool adaptolbetter_;  ///< adaptive tolerance for solver useful?
    Teuchos::RCP<Core::LinAlg::Solver> solver_;  ///< solver for linear standard linear system
    int counter_;                                ///< counts how often #Solve is called
    Inpar::Solid::ConSolveAlgo algochoice_;
    Inpar::Solid::StcScale stcalgo_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat_;

  };  // class
}  // namespace CONSTRAINTS
FOUR_C_NAMESPACE_CLOSE

#endif
