/*!----------------------------------------------------------------------
\file uzawasolver.H

\brief Class constaining uzawa algorithm to solve linear system.

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef UZAWASOLVER_H
#define UZAWASOLVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>


#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"


using namespace std;
using namespace Teuchos;

namespace UTILS
{
  /*!
  \brief Class containing uzawa algorithm to solve linear system.
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class UzawaSolver
  {
  public:

    //! @name Constructors and destructors and related methods

    /*!
    \brief Constructor
    */

    UzawaSolver(
        RCP<DRT::Discretization> discr,     ///< Discretization
        LINALG::Solver& solver,             ///< Solver to solve linear subproblem in iteration
        RCP<Epetra_Vector>    dirichtoggle, ///< Vector used for ApplyDirichletToSystem
        RCP<Epetra_Vector>    invtoggle,    ///< Vector to blank residuals of Dirichlet DOFs
        ParameterList param);               ///< parameterlist containing solver parameters

    /*!
        \brief Destructor

     */
    ~UzawaSolver()
    {
        return;
    };

    //@}

    //! @name Public Member Functions
    /*!
      \brief Solve linear system using uzawa algorithm to deal with zero entries on the diagonal;

    */
    void Solve(
        RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
        RCP<LINALG::SparseMatrix> constr, ///< constraint matrix
        RCP<Epetra_Vector> disp_inc,      ///< displacement increment to compute
        RCP<Epetra_Vector> disp_lagr,     ///< lagrange multiplier increment to compute
        RCP<Epetra_Vector> rhs_standard,  ///< standard rhs
        RCP<Epetra_Vector> rhs_constr     ///< constraint errors
        );

    /*!
       \brief Return the current value of the uzawa parameter
    */
    double GetUzawaParameter()
    {
      return uzawaParam_;
    }

    /*!
      \brief Set uzawa parameter (used for dynamic adaptation and restart)
    */
    void SetUzawaParameter(double restartval)
    {
      uzawaParam_=restartval;
      counter_=1;
      return;
    }

    //@}

private:

    // do not want = operator, cctor

    UzawaSolver operator = (const UzawaSolver& old);
    UzawaSolver(const UzawaSolver& old);


    //! name@ Private Attributes
    RCP<DRT::Discretization> actdisc_;
    int maxIter_;
    int tolerance_;
    bool tolcheck_;
    double uzawaParam_;
    double norm_uzawa_;
    double norm_uzawa_alt_;
    double quotient_;
    double norm_constr_uzawa_;
    double toldisp_;
    double tolres_;
    double tolconstr_;
    RCP<Epetra_Vector>    dirichtoggle_;
    RCP<Epetra_Vector>    invtoggle_;
    bool isadapttol_;
    bool adaptolbetter_ ;
    RCP<LINALG::Solver> solver_;
    int counter_;
    //@}


  }; //class
}
#endif /*CONSTRAINTMANAGER_H*/
