/*!----------------------------------------------------------------------
\file avms_solver.H

\class VM3_Solver

\brief An algebraic variational multiscale solver and method

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef AVMS_SOLVER_H
#define AVMS_SOLVER_H

#ifdef CCADISCRET

// Trilinos includes
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_IntVector.h"
#include "Epetra_CrsMatrix.h"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_genP.h"
#include "ml_operator.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"


using namespace std;
using namespace Teuchos;
using namespace MLAPI;


/*!
\class AVMS_Solver

\brief An algebraic variational multiscale solver and method

The AVMS (formerly VM3) method class. It takes 2 matrices upon construction:
Aforfine is a fine resolved scale system that includes stabilization to 
be present on the fine resolved scales.
Aforcoarse is a fine resolved scale system that includes stabilization that
is supposd to be present on the algebraic coarse resolved scale system as well.
This class will build a coarse resolved scale equation from Aforcoarse using
plain aggregation an will create (and later solve) a block 2x2 system of fine and coarse
resolved scales:

\code

    -                                    -   -      -      -          -
   |                                     |  |       |     |           |
   | P^T Aforcoarse P    P^T Aforcoarse  |  | x_c   |     |  P^T f_f  |
   |                                     |  |       |  =  |           |
   |  Aforcoarse P         Aforfine      |  | dx_f  |     |  f_f      |
   |                                     |  |       |     |           |
    -                                    -   -      -      -          -

\endcode

The solution is returned as

\code

   x = dx_f + P x_c

\endcode

\note At no point, the user needs to supply coarse level quantities.

\author vg/gee

*/

class  AVMS_Solver
{
public:

  /*!
  \brief ctor

  \param Aforfine (in)  : fine resolved scale equation includes fine scale stabilization
  \param Aforcoarse (in): fine resolved scale equation the algebraic coarse scale is build from
                          (which exludes some of the stabilization present in Aforfine)
  \param mlparams (in/out): valid ML parameter list, the nullspace is altered by this class
  */
  explicit AVMS_Solver(RCP<Epetra_CrsMatrix> Aforfine,
                       RCP<Epetra_CrsMatrix> Aforcoarse,
                       ParameterList& mlparams);

  /*!
  \brief Destroys an instance of this class

  Destructor

  */
  virtual ~AVMS_Solver();

  /*!
  \brief Apply the avms solver

  */
  int Solve(const Epetra_Vector& X, Epetra_Vector& Y, ParameterList& params);



private:
  // don't want copy-ctor and = operator
  AVMS_Solver(AVMS_Solver& old);
  AVMS_Solver operator = (const AVMS_Solver& old);

  /*!
  \brief Compute the avms setup

  */
  bool Compute(RCP<Epetra_CrsMatrix> Aforfine,
               RCP<Epetra_CrsMatrix> Aforcoarse);

  ParameterList&         mlparams_;  // ML parameter list
  RCP<Epetra_CrsMatrix>  Acombined_; // input matrix without fine-scale term
  RCP<Epetra_CrsMatrix>  P_;         // prolongation from algerbaic coarse level to this level

}; // AVMS_Solver


#endif

#endif
