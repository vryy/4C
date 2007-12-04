#ifndef VM3_SOLVER_H
#define VM3_SOLVER_H

#ifdef CCADISCRET

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

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
\class VM3_Solver

\brief


\author vg

*/

class  VM3_Solver
{
public:

  // @{ \name Constructors and destructors


  /*!
  \brief ctor

  */
  explicit VM3_Solver( RefCountPtr<Epetra_CrsMatrix>& Aplus,
                       RefCountPtr<Epetra_CrsMatrix>& A,
                       RefCountPtr<Epetra_Vector>& rplus,
                       RefCountPtr<Epetra_Vector>& r,
                       RefCountPtr<Epetra_Vector>& sol,
                       const RefCountPtr<Epetra_Vector> dbctoggle,
                       ParameterList& mlparams,
                       bool increm);

  /*!
  \brief Destroys an instance of this class

  Destructor

  */
  virtual ~VM3_Solver() {}

  //@}
  // @{ \name Query methods

  /*!
  \brief get Comm of this class

  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

  */
  const Epetra_Comm & Comm() const {return(A_->Comm());}


  //@}
  // @{ \name Apply methods


  /*!
  \brief Apply the VM3 solver

   Perform one sweep of the VM3 solver


   \param X  (In)  : Epetra_MultiVector matching the fine level map of this
                     solver
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int Solve(const Epetra_Vector& X, Epetra_Vector& Y, ParameterList& params);


  /*!
  \brief Apply the VM3 solver for incremental formulation

   Perform one sweep of the VM3 solver for incremental formulation


   \param X  (In)  : Epetra_MultiVector matching the fine level map of this
                     solver
   \param X2  (In)  : modified Epetra_MultiVector matching the fine level map of this
                     solver
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int Solve(const Epetra_Vector& X, const Epetra_Vector& X2, Epetra_Vector& Y, ParameterList& params);


private:

  // don't want copy-ctor and = operator
  VM3_Solver(VM3_Solver& old);
  VM3_Solver operator = (const VM3_Solver& old);

  bool                           increm_;  // flag indicating incremental formulation 

  ParameterList&                 mlparams_;    // ML parameter list

  RCP<Epetra_CrsMatrix>& Aplus_;     // input matrix including fine-scale term
  RCP<Epetra_CrsMatrix>& A_;         // input matrix without fine-scale term
  RCP<Epetra_Vector>& rplus_;
  RCP<Epetra_Vector>& r_;
  RCP<Epetra_Vector>& sol_;
  //RCP<Epetra_CrsMatrix>  Acombined_; // input matrix without fine-scale term
  Operator               Ptent_;
  Operator               Rtent_;
  //RCP<Epetra_Map>        coarsermap_;

}; // VM3_Solver


#endif

#endif
