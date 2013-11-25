/*!----------------------------------------------------------------------
\file windkesselsolver.H

\brief Class containing direct solver to solve linear system.

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#ifndef WINDKESSELSOLVER_H
#define WINDKESSELSOLVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_structure/strtimint.H"

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
  class WindkesselManager;
}


namespace UTILS
{
  /*!
  \brief Class containing direct solver to solve linear system.

  \author Hirschvogel (hirschvogel@lnm.mw.tum.de)
  */
  class WindkesselSolver
  {
  public:

    /*!
    \brief Constructor
    */
    WindkesselSolver
    (
      RCP<DRT::Discretization> discr,     ///< Discretization
      LINALG::Solver& solver,             ///< Solver to solve linear subproblem in iteration
      RCP<LINALG::MapExtractor> dbcmaps,  ///< Map extractor for Dirichlet DOFs
      ParameterList param                 ///< parameterlist containing solver parameters
    );

    /*!
        \brief Destructor
     */
    ~WindkesselSolver()
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
      \brief Solve Windkessel linear system
    */
    void Solve
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> windkoffdiag1, ///< Windkessel matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> windkoffdiag2, ///< transpose of Windkessel matrix without Dirichlet zeros
      RCP<LINALG::SparseMatrix> windkstiff, ///< Windkessel stiffness
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> disppres,     ///< pressure increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhswindk     ///< Windkessel errors
    );

    void SetSTCProp
    (
      INPAR::STR::STC_Scale stcalgo,
      RCP<LINALG::SparseMatrix> stcmat
    )
    {
      stcalgo_=stcalgo;
      stcmat_=stcmat;
    };

    //! Return Teuchos::rcp to WindkesselManager windkman_
    const RCP< ::UTILS::WindkesselManager> GetWindkesselManager()
    {
      return windkman_;
    }

private:

    // do not want = operator, cctor
    WindkesselSolver operator = (const WindkesselSolver& old);
    WindkesselSolver(const WindkesselSolver& old);
    Teuchos::RCP<UTILS::WindkesselManager> windkman_;  //!< Windkessel manager
    
    /*!
      \brief Solve linear system directly by assembling everything into one big matrix
    */
    void SolveDirect
    (
      RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      RCP<LINALG::SparseMatrix> windkoffdiag1, ///< Windkessel matrix with Dirichlet zeros
      RCP<LINALG::SparseMatrix> windkoffdiag2, ///< transpose of Windkessel matrix without Dirichlet zeros
      RCP<LINALG::SparseMatrix> windkstiff, ///< Windkessel stiffness
      RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      RCP<Epetra_Vector> disppres,     ///< pressure increment to compute
      const RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const RCP<Epetra_Vector> rhswindk     ///< Windkessel errors
    );

    RCP<DRT::Discretization> actdisc_; ///<standard discretization
    int maxIter_;  ///< number of maximal iterations
    double iterationparam_; ///< parameter for Uzawa algorithm
    double minparam_; ///< minimal possible parameter for Uzawa algorithm
    double iterationtol_;  ///< tolerance
    double tolres_;  ///< tolerance for residual
    double tolwindk_;  ///< tolerance for Windkessel
    RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    RCP<Epetra_Vector>    firstdispinc_; ///< history variable holding displacement increment for first NRI
    RCP<Epetra_Vector>    firstpresinc_; ///< history variable holding pressure increment for first NRI
    bool isadapttol_;  ///< adaptive tolerance for solver?
    bool adaptolbetter_ ; ///< adaptive tolerance for solver useful?
    RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    INPAR::STR::ConSolveAlgo algochoice_;
    INPAR::STR::STC_Scale stcalgo_;
    RCP<LINALG::SparseMatrix> stcmat_;

  }; //class
}
#endif /*WINDKESSELSOLVER_H*/
