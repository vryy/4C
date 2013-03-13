/*!----------------------------------------------------------------------
\file constraintsolver.H

\brief Class containing direct solver and uzawa algorithm to solve linear system.

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef CONSTRAINTSOLVER_H
#define CONSTRAINTSOLVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_structure.H"

using namespace Teuchos;

//forward declarations
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
}

namespace UTILS
{
  /*!
  \brief Class containing uzawa algorithm to solve linear system.

  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class ConstraintSolver
  {
  public:

    /*!
    \brief Constructor
    */
    ConstraintSolver
    (
      RCP<DRT::Discretization> discr,     ///< Discretization
      LINALG::Solver& solver,             ///< Solver to solve linear subproblem in iteration
      RCP<LINALG::MapExtractor> dbcmaps,  ///< Map extractor for Dirichlet DOFs
      ParameterList param                 ///< parameterlist containing solver parameters
    );

    /*!
        \brief Destructor
     */
    ~ConstraintSolver()
    {
        return;
    };

    /*!
      \brief Set it up
     */
    void Setup
    (
      RCP<DRT::Discretization> discr,
      LINALG::Solver& solver,
      RCP<LINALG::MapExtractor> dbcmaps,
      ParameterList params
    );

    /*!
      \brief Solve constraint linear system
    */
    void Solve
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> constr, ///< constraint matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> constrT, ///< transpose of constraint matrix without Dirichlet zeros
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> displagr,     ///< lagrange multiplier increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    /*!
       \brief Return the current value of the uzawa parameter
    */
    double GetUzawaParameter()
    {
      return iterationparam_;
    }

    /*!
      \brief Set uzawa parameter (used for dynamic adaptation and restart)
    */
    void SetUzawaParameter
    (
      double restartval ///< value to replace Uzawa Parameter with
    )
    {
      iterationparam_=restartval;
      counter_=1;
      return;
    }

    void SetSTCProp
    (
      INPAR::STR::STC_Scale stcalgo,
      RCP<LINALG::SparseMatrix> stcmat
    )
    {
      stcalgo_=stcalgo;
      stcmat_=stcmat;
    };

private:

    // do not want = operator, cctor
    ConstraintSolver operator = (const ConstraintSolver& old);
    ConstraintSolver(const ConstraintSolver& old);

    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;
    */
    void SolveUzawa
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> constr, ///< constraint matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> constrT, ///< transpose of constraint matrix without Dirichlet zeros
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> displagr,     ///< lagrange multiplier increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );
    
    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;
    */
    void SolveSimple
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> constr, ///< constraint matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> constrT, ///< transpose of constraint matrix without Dirichlet zeros
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> displagr,     ///< lagrange multiplier increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );
    
    /*!
      \brief Solve linear system directly by assembling everything into one big matrix
    */
    void SolveDirect
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> constr, ///< constraint matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> constrT, ///< transpose of constraint matrix without Dirichlet zeros
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> displagr,     ///< lagrange multiplier increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhsconstr     ///< constraint errors
    );

    RCP<DRT::Discretization> actdisc_; ///<standard discretization
    int maxIter_;  ///< number of maximal iterations
    double iterationparam_; ///< parameter for Uzawa algorithm
    double minparam_; ///< minimal possible parameter for Uzawa algorithm
    double iterationtol_;  ///< tolerance
    double tolres_;  ///< tolerance for residual
    double tolconstr_;  ///< tolerance for constraint
    RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    RCP<Epetra_Vector>    firstdispinc_; ///< history variable holding displacement increment for first NRI
    RCP<Epetra_Vector>    firstlagrinc_; ///< history variable holding multiplier increment for first NRI
    bool isadapttol_;  ///< adaptive tolerance for solver?
    bool adaptolbetter_ ; ///< adaptive tolerance for solver useful?
    RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    INPAR::STR::ConSolveAlgo algochoice_;
    INPAR::STR::STC_Scale stcalgo_;
    RCP<LINALG::SparseMatrix> stcmat_;

  }; //class
}
#endif /*CONSTRAINTMANAGER_H*/
