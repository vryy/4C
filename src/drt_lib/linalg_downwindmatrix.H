/*!----------------------------------------------------------------------
\file linalg_downwindmatrix.H

\class LINALG::DownwindMatrix


<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_DOWNWINDMATRIX_H
#define LINALG_DOWNWINDMATRIX_H

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_systemmatrix.H"
#include "linalg_ana.H"

using namespace std;
using namespace Teuchos;

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class DownwindMatrix


\author gee (gee@lnm.mw.tum.de)
*/
class  DownwindMatrix : public virtual Epetra_RowMatrix
{
public:

  /*!
  \brief Standard Constructor
  */
  explicit DownwindMatrix(RCP<Epetra_CrsMatrix> A);

  /*!
  \brief Destructor
  */
  virtual ~DownwindMatrix() {}

  /*!
  \brief Returns the label of this class.
  */
  const char* Label() const { return "LINALG::DownwindMatrix"; }

  /*!
  \brief get Comm of this class

  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

  */
  const Epetra_Comm& Comm() const {return(A_->Comm());}


  /*!
  \brief Get fine level OperatorDomainMap

  Derived from Epetra_Operator, get fine level OperatorDomainMap

  */
  const Epetra_Map& OperatorDomainMap() const {return A_->OperatorDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap

  Derived from Epetra_Operator, get fine level OperatorRangeMap

  */
  const Epetra_Map& OperatorRangeMap() const {return A_->OperatorRangeMap();}

  /*!
  \brief ApplyInverse the preconditioner

   ApplyInverse the preconditioner. Method is derived from Epetra_Operator.


   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  { return A_->ApplyInverse(X,Y); }

  /*!
  \brief not implemented
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  { return A_->Apply(X,Y); }

  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose)
  { return A_->SetUseTranspose(UseTranspose); }

  /*!
  \brief not implemented
  */
  double NormInf() const
  { return A_->NormInf(); }


  /*!
  \brief not implemented
  */
  bool UseTranspose() const
  { return A_->UseTranspose(); }

  /*!
  \brief not implemented
  */
  bool HasNormInf() const
  { return A_->HasNormInf(); }

  //! Returns a reference to the Epetra_BlockMap for this object.
  virtual const Epetra_BlockMap & Map() const
  { return A_->Map(); }

    //! Returns the number of nonzero entries in MyRow.
    /*! 
    \param In
           MyRow - Local row.
    \param Out
	   NumEntries - Number of nonzero values present.
	  
    \return Integer error code, set to 0 if successful.
    */
    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
    {  return A_->NumMyRowEntries(MyRow,NumEntries); }

    //! Returns the maximum of NumMyRowEntries() over all rows.
    virtual int MaxNumEntries() const
    {  return A_->MaxNumEntries(); }
    
    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           MyRow - Local row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
    { return A_->ExtractMyRowCopy(MyRow,Length,NumEntries,Values,Indices); }
    
    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
    { return A_->ExtractDiagonalCopy(Diagonal); }

    //! Returns the result of a Epetra_RowMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { return A_->Multiply(TransA,X,Y); }

    //! Returns result of a local-only solve using a triangular Epetra_RowMatrix with Epetra_MultiVectors X and Y.
    /*! This method will perform a triangular solve independently on each processor of the parallel machine.
        No communication is performed.
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
                      Epetra_MultiVector& Y) const
    { dserror("LINALG::DownwindMatrix does not implement Epetra_RowMatrix::Solve"); return -1; }

    //! Computes the sum of absolute values of the rows of the Epetra_RowMatrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
    */
    virtual int InvRowSums(Epetra_Vector& x) const
    { return A_->InvRowSums(x); }

    //! Computes the sum of absolute values of the columns of the Epetra_RowMatrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
    */
    virtual int InvColSums(Epetra_Vector& x) const
    { return A_->InvColSums(x); }

    //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
    */
    virtual int RightScale(const Epetra_Vector& x)
    { return A_->RightScale(x); }

    //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
    */
    virtual int LeftScale(const Epetra_Vector& x)
    { return A_->LeftScale(x); }

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    virtual bool Filled() const
    { return A_->Filled(); }

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    virtual double NormOne() const
    { return A_->NormOne(); }

    //! Returns the number of nonzero entries in the global matrix.
    /*
      Note that depending on the matrix implementation, it is sometimes
      possible to have some nonzeros that appear on multiple processors.
      In that case, those nonzeros may be counted multiple times (also
      depending on the matrix implementation).
    */
    virtual int NumGlobalNonzeros() const
    { return A_->NumGlobalNonzeros(); }

    //! Returns the number of global matrix rows.
    virtual int NumGlobalRows() const
    { return A_->NumGlobalRows(); }

    //! Returns the number of global matrix columns.
    virtual int NumGlobalCols() const
    { return A_->NumGlobalCols(); }

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumGlobalDiagonals() const
    { return A_->NumGlobalDiagonals(); }
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    virtual int NumMyNonzeros() const
    { return A_->NumMyNonzeros(); }

    //! Returns the number of matrix rows owned by the calling processor.
    virtual int NumMyRows() const
    { return A_->NumMyRows(); }

    //! Returns the number of matrix columns owned by the calling processor.
    virtual int NumMyCols() const
    { return A_->NumMyCols(); }

    //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
    virtual int NumMyDiagonals() const
    { return A_->NumMyDiagonals(); }

    //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool LowerTriangular() const
    { return A_->LowerTriangular(); }

    //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
    virtual bool UpperTriangular() const
    { return A_->UpperTriangular(); }

    //! Returns the Epetra_Map object associated with the rows of this matrix.
    virtual const Epetra_Map & RowMatrixRowMap() const
    { return A_->RowMatrixRowMap(); }

    //! Returns the Epetra_Map object associated with the columns of this matrix.
    virtual const Epetra_Map & RowMatrixColMap() const
    { return A_->RowMatrixColMap(); }

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const
    { return A_->RowMatrixImporter(); }




private:

  // don't want copy-ctor and = operator
  DownwindMatrix(DownwindMatrix& old);
  DownwindMatrix operator = (const DownwindMatrix& old);

  /*!
  \brief setup phase of preconditioner
  */
  void Setup(RCP<Epetra_Operator> A);


  RCP<Epetra_CrsMatrix> A_;

}; // class  DownwindMatrix : public virtual Epetra_Operator
}  // namespace LINALG


#endif // LINALG_DOWNWINDMATRIX_H
#endif  // #ifdef CCADISCRET
