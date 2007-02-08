/*!----------------------------------------------------------------------
\file linalg_solver.H

\class LINALG::Solver

\brief A general solver interface to Trilinos solvers and spooles

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef LINALG_SOLVER_H
#define LINALG_SOLVER_H

#include "vector"
#include "string"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;


// forward declarations
struct _SOLVAR;


/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class Solver

\brief A general solver interface to Trilinos solvers and spooles

\author gee (gee@lnm.mw.tum.de)
*/
class Solver
{
public:


  /*!
  \brief Standard Constructor
  
  */
  Solver();

  /*!
  \brief Copy Constructor

  */
  Solver(const LINALG::Solver& old);

  /*!
  \brief Destructor

  */
  virtual ~Solver();

  /*!
  \brief Translate solver input parameters from old ccarat style to
         new solver parameters list style
         
  \params params (out): parameter list holding new style solver parameters
                        on output
  \params actsolv (in): ccarat style struct holding solver parameters

  */
  void TranslateSolverParameters(ParameterList& params, struct _SOLVAR* actsolv);

protected:

private:

  // don't want = operator
  Solver operator = (const Solver& old);


}; // class Solver
} // namespace LINALG











#endif  // #ifndef LINALG_SOLVER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
