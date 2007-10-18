/*!
 * \file mfoperator.H
 *
 * \class MatrixFreeOperator
 *
 * \brief Approximation to matrix-vector product
 *
 * \date 10/07
 *
 */

#ifdef CCADISCRET

#ifndef MFOPERATOR_H
#define MFOPERATOR_H

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"   // for testing only!!!!

// Teuchos
#include "Teuchos_RefCountPtr.hpp"

//#include "strugenalpha.H"

// ---------- Forward Declarations ----------
using namespace Teuchos;
class StruGenAlpha;

namespace LINALG
{

/*!
\brief First-Order Taylor series expansion approximation to
matrix-vector product

\author  L. Wiechert (wiechert@lnm.mw.tum.de)
*/

class MatrixFreeOperator : public virtual Epetra_Operator
{
public:

  explicit MatrixFreeOperator(StruGenAlpha& integrator, RCP<Epetra_CrsMatrix> stiff);

  //! Destroys the operator.
  /*!
   Destroys the MatrixFreeOperator.
   */
  virtual ~MatrixFreeOperator() {}

  //! Returns the label of this class.
  /*!
   Returns the label of this class.
   */
  const char* Label() const {return(&(label_[0]));}

  //! get Comm of this class
  /*!
   Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class
   */
  const Epetra_Comm & Comm() const {return(du_.Comm());}

  //! Get OperatorDomainMap.
  /*!
   Derived from Epetra_Operator, get fine level OperatorDomainMap.
   */
  const Epetra_Map & OperatorDomainMap() const;

  //! Get OperatorRangeMap.
  /*!
   Derived from Epetra_Operator, get fine level OperatorRangeMap.
   */
  const Epetra_Map & OperatorRangeMap() const;

  //! ApplyInverse the operator, not implemented
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  { cout << "MatrixFreeOperator::ApplyInverse does not make sense";
    return(-1);}

  //! Apply the operator
  /*!
   Approximate calculation of the matrix-vector product (first-order
   Taylor series expansion).
   */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! SetUseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  int SetUseTranspose(bool UseTranspose) { return -1; }

  //! NormInf, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  double NormInf() const {return(-1.0);}

  //! UseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  bool UseTranspose() const {return false;}

  //! HasNormInf, not implemented.
  /*!
   Derived from Epetra_Operator, always returns false.
   */
  bool HasNormInf() const {return(false);}

private:

  explicit MatrixFreeOperator(MatrixFreeOperator& old);
  MatrixFreeOperator operator = (const MatrixFreeOperator& old);

  string                                                    label_;
  StruGenAlpha&                                        integrator_;
  Epetra_Vector                                                du_;
  Epetra_Vector                                                 F_;
  RCP<Epetra_CrsMatrix>                                     stiff_;
  const Epetra_Map*                                           map_;

}; // class JFNK_Operator
}  // namespace LINALG

#endif
#endif
