/*!----------------------------------------------------------------------
\file linalg_mlapi_operator.H

\class LINALG::AMG_Operator

\brief A multipurpose experimental multigrid operator

This operator based on the ml advanced programming interface is a
multipurpose development object for amg ideas that shall be tested in
the baci framework

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef LINALG_MLAPI_OPERATOR_H
#define LINALG_MLAPI_OPERATOR_H

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_genP.h"
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

#include "drt_dserror.H"

using namespace std;
using namespace Teuchos;
using namespace MLAPI;

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{
/*!
\class AMG_Operator

\brief A multipurpose experimental multigrid operator

This operator based on the ml advanced programming interface is a
multipurpose development object for amg ideas that shall be tested in
the baci framework

\author gee (gee@lnm.mw.tum.de)
*/
class  AMG_Operator : public virtual Epetra_Operator
{
public:

  /*!
  \brief Standard Constructor
  */
  explicit AMG_Operator(RCP<Epetra_CrsMatrix> A, 
                        ParameterList& params, 
                        const bool compute);

  /*!
  \brief Destructor
  */
  virtual ~AMG_Operator() {}

  /*!
  \brief Returns the label of this class.
  */
  const char* Label() const;// {return(&(label_[0]));}

  /*!
  \brief get Comm of this class
  
  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class 
  
  */
  const Epetra_Comm & Comm() const {return(A_->Comm());}

  /*!
  \brief Get fine level OperatorDomainMap
  
  Derived from Epetra_Operator, get fine level OperatorDomainMap
  
  */
  const Epetra_Map & OperatorDomainMap() const {return A_->OperatorDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap
  
  Derived from Epetra_Operator, get fine level OperatorRangeMap
  
  */
  const Epetra_Map & OperatorRangeMap() const {return A_->OperatorRangeMap();}
  
  /*!
  \brief ApplyInverse the preconditioner
  
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator. 
   
   
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  /*!
  \brief not implemented
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { 
    dserror("Apply does not make sense for LINALG::AMG_Operator");
    return(-1);
  }
    
  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose) 
  { 
    dserror("SetUseTranspose not impl.");
    return -1;
  }

  /*!
  \brief not implemented
  */
  double NormInf() const 
  {
    dserror("NormInf not impl.");
    return(-1.0);
  }


  /*!
  \brief not implemented
  */
  bool UseTranspose() const 
  {
    dserror("UseTranspose not impl.");
    return false;
  }

  /*!
  \brief not implemented
  */
  bool HasNormInf() const 
  {
    dserror("HasNormInf not impl.");
    return false;
  }

private:

  /*!
  \brief Compute the preconditioner
  
   Compute the mortar modified ml preconditioner
  */
  void Compute();  

  // don't want copy-ctor and = operator
  AMG_Operator(AMG_Operator& old);
  AMG_Operator operator = (const AMG_Operator& old);
  
  string                label_;  // name of this operator
  ParameterList&        params_; // ml parameter list (intput)
  RCP<Epetra_CrsMatrix> A_;      // input matrix  



}; // class  AMG_Operator : public virtual Epetra_Operator
}  // namespace LINALG


#endif // LINALG_MLAPI_OPERATOR_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
