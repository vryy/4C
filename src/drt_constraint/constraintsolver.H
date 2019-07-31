/*!----------------------------------------------------------------------
\brief Class containing Uzawa algorithm to solve linear system.
\level 2
\maintainer Matthias Mayr
*----------------------------------------------------------------------*/

#ifndef CONSTRAINTSOLVER_H
#define CONSTRAINTSOLVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MapExtractor;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace UTILS
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
    ConstraintSolver(Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization
        LINALG::Solver& solver,  ///< Solver to solve linear subproblem in iteration
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps,  ///< Map extractor for Dirichlet DOFs
        Teuchos::ParameterList param                 ///< parameterlist containing solver parameters
    );

    /*!
        \brief Destructor
     */
    ~ConstraintSolver() { return; };

    /*!
      \brief Set it up
     */
    void Setup(Teuchos::RCP<DRT::Discretization> discr, LINALG::Solver& solver,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, Teuchos::ParameterList params);

    /*!
      \brief Solve constraint linear system
    */
    void Solve(Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
        Teuchos::RCP<LINALG::SparseMatrix> constr,  ///< constraint matrix with Dirichlet zeros
        Teuchos::RCP<LINALG::SparseMatrix>
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Epetra_Vector> dispinc,   ///< displacement increment to compute
        Teuchos::RCP<Epetra_Vector> displagr,  ///< lagrange multiplier increment to compute
        const Teuchos::RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
        const Teuchos::RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    /*!
       \brief Return the current value of the uzawa parameter
    */
    double GetUzawaParameter() { return iterationparam_; }

    /*!
      \brief Set uzawa parameter (used for dynamic adaptation and restart)
    */
    void SetUzawaParameter(double restartval  ///< value to replace Uzawa Parameter with
    )
    {
      iterationparam_ = restartval;
      counter_ = 1;
      return;
    }

    void SetSTCProp(INPAR::STR::STC_Scale stcalgo, Teuchos::RCP<LINALG::SparseMatrix> stcmat)
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
    void SolveUzawa(Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
        Teuchos::RCP<LINALG::SparseMatrix> constr,  ///< constraint matrix with Dirichlet zeros
        Teuchos::RCP<LINALG::SparseMatrix>
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Epetra_Vector> dispinc,   ///< displacement increment to compute
        Teuchos::RCP<Epetra_Vector> displagr,  ///< lagrange multiplier increment to compute
        const Teuchos::RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
        const Teuchos::RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;
    */
    void SolveSimple(Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
        Teuchos::RCP<LINALG::SparseMatrix> constr,  ///< constraint matrix with Dirichlet zeros
        Teuchos::RCP<LINALG::SparseMatrix>
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Epetra_Vector> dispinc,   ///< displacement increment to compute
        Teuchos::RCP<Epetra_Vector> displagr,  ///< lagrange multiplier increment to compute
        const Teuchos::RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
        const Teuchos::RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    /*!
      \brief Solve linear system directly by assembling everything into one big matrix
    */
    void SolveDirect(Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
        Teuchos::RCP<LINALG::SparseMatrix> constr,  ///< constraint matrix with Dirichlet zeros
        Teuchos::RCP<LINALG::SparseMatrix>
            constrT,  ///< transpose of constraint matrix without Dirichlet zeros
        Teuchos::RCP<Epetra_Vector> dispinc,   ///< displacement increment to compute
        Teuchos::RCP<Epetra_Vector> displagr,  ///< lagrange multiplier increment to compute
        const Teuchos::RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
        const Teuchos::RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    Teuchos::RCP<DRT::Discretization> actdisc_;  ///< standard discretization
    int maxIter_;                                ///< number of maximal iterations
    double iterationparam_;                      ///< parameter for Uzawa algorithm
    double minparam_;                            ///< minimal possible parameter for Uzawa algorithm
    double iterationtol_;                        ///< tolerance
    double tolres_;                              ///< tolerance for residual
    double tolconstr_;                           ///< tolerance for constraint
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;   ///< \b only for compatability: dirichlet toggle --
                                                 ///< monitor its target change!
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Epetra_Vector>
        firstdispinc_;  ///< history variable holding displacement increment for first NRI
    Teuchos::RCP<Epetra_Vector>
        firstlagrinc_;    ///< history variable holding multiplier increment for first NRI
    bool isadapttol_;     ///< adaptive tolerance for solver?
    bool adaptolbetter_;  ///< adaptive tolerance for solver useful?
    Teuchos::RCP<LINALG::Solver> solver_;  ///< solver for linear standard linear system
    int counter_;                          ///< counts how often #Solve is called
    INPAR::STR::ConSolveAlgo algochoice_;
    INPAR::STR::STC_Scale stcalgo_;
    Teuchos::RCP<LINALG::SparseMatrix> stcmat_;

  };  // class
}  // namespace UTILS
#endif /*CONSTRAINTMANAGER_H*/
