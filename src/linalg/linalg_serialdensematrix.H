/*!----------------------------------------------------------------------
\file linalg_serialdensematrix.H

\brief Declaration of serial dense matrix wrapper class

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef SERIALDENSEMATRIX_H
#define SERIALDENSEMATRIX_H


#include "Epetra_SerialDenseMatrix.h"


namespace LINALG
{

/*!
\brief A class that wraps Epetra_SerialDenseMatrix with minor modifications
       in the constructor and a modified memory management

Epetra_SerialDenseMatrix always allocates new memory upon call of Shape() or
Reshape(). This class re-uses previously allocated memory as much as possible,
similar to the implementation in std::vector. If the matrix size is reduced
or if the allocated memory is sufficient for the new matrix sizes, the old
data array is re-used. Re-allocation only happens in case more memory is
necessary. This way, some memory might be wasted but it results in code that
is generally considerably faster and avoids memory fragmentation. To force a
new allocation with the correct amount of memory, call Shape(0,0) first.

*/
class SerialDenseMatrix : public Epetra_SerialDenseMatrix
{
public:

  //! @name Constructors and destructors and related methods
  //@{

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseMatrix(bool set_object_label=true);

  */
  SerialDenseMatrix(bool set_object_label=true);

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseMatrix(int NumRows, int NumCols, bool set_object_label=true);

  Opposed to the base class constructor, this one does NOT initialize the matrix to zero.

  \param NumRows          (in): no. of rows the matrix is supposed to have
  \param NumCols          (in): no. of cols the matrix is supposed to have
  \param init             (in): flag indicating whether matrix is zeroed out (defaults to false!)
  \param set_object_label (in): set an epetra object label for print out (defaults to true!)

  */
  SerialDenseMatrix(int NumRows, int NumCols, bool init=false, bool set_object_label=true);

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA, int NumRows, int NumCols,
                                  bool set_object_label=true);

  \param CV               (in): specifies whether a view to the data in A is created or it is copied
  \param A                (in): pointer to data array, organized by columns
  \param LDA              (in): leading dimension of A, i.e., the distance between columns in the data array.
                                The first column starts at A, the second at A+LDA, the third at A+2*LDA, etc.
  \param NumRows          (in): no. of rows the matrix is supposed to have, must be smaller or equal LDA
  \param NumCols          (in): no. of cols the matrix is supposed to have
  \param set_object_label (in): set an epetra object label for print out (defaults to true!)

  */
  SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA, int NumRows, int NumCols,
                           bool set_object_label=true);

  /*!
  \brief Standard Constructor allowing view on Epetra_SerialDenseMatrix and wrapping
         SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA, int NumRows, int NumCols,
                           bool set_object_label=true);


  */
  SerialDenseMatrix(Epetra_SerialDenseMatrix& Source, Epetra_DataAccess CV, bool set_object_label=true);

  /*!
  \brief Standard Copy Constructor wraps
         Epetra_SerialDenseMatrix(const SerialDenseMatrix& Source);

  */
  SerialDenseMatrix(const SerialDenseMatrix& Source);

  /*!
  \brief Standard Copy Constructor wraps
         Epetra_SerialDenseMatrix(const SerialDenseMatrix& Source);

  */
  SerialDenseMatrix(const Epetra_SerialDenseMatrix& Source);

  /*!
  \brief Shape dimensions of this matrix as in Epetra_SerialDenseMatrix::Shape
         (but only allocate new memory if necessary)

  If the number of rows or columns is zero, the memory will be freed. Otherwise,
  the algorithm keeps the old memory if it is sufficient and only sets the
  matrix bounds differently.
  */
  int Shape(int NumRows, int NumCols);

  /*!
  \brief Shape dimensions of this matrix but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseMatrix::Shape(...))
  */
  int LightShape(int NumRows, int NumCols);

  /*!
  \brief Reshape dimensions of this matrix and zero out excess space

  Existing values are kept as far as they fit into the new shape and
  potential excess space in the matrix is zeroed out

  If the number of rows or columns is zero, the memory will be freed. Otherwise,
  the algorithm keeps the old memory if it is sufficient and only sets the
  matrix bounds differently.
  */
  int Reshape(int NumRows, int NumCols);

  /*!
  \brief Reshape dimensions of this matrix but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseMatrix::Reshape(...))

  Existing values are kept as far as they fit into the new shape but
  potential excess space in the matrix is NOT zeroed out

  */
  int LightReshape(int NumRows, int NumCols);

  /*!
  \brief Destructor

  */
  virtual ~SerialDenseMatrix();

  /*!
  \brief Assignment operator from another SerialDenseMatrix

  If the matrix memory is enough to keep the new matrix, no memory reallocation
  is done, similarly to memory management of the LightShape method.
   */
  SerialDenseMatrix & operator= (const SerialDenseMatrix &Source);

  //@}

  //! @name Mathematical methods
  //@{

  /*!
  \brief Determinant Computation using the Sarrus rule.
  Internal computation is based on the long double data type (80/128 bit depending on platform)
  this allows for higher precision when used on a bigger matrix. Long double is also used for output,
  therefore it has to be expicitly casted to a double if used in a "double only" context.

  \author volf (alexander.volf@mytum.de)
  */
  long double Det_long();

  /*!
   * \brief Update matrix components with scaled values of B,
   *        this = ScalarThis * this + ScalarB * B
   * \author bborn \date 08/08
   */
  void Update(
    const double& ScalarB,  /*!< scale input matrix with */
    const Epetra_SerialDenseMatrix& B,  /*!< input matrix */
    const double& ScalarThis  /*!< scale this matrix with */
  );

  /*!
   * \brief Set matrix components to zero,
   *        this = 0.0
   * \author bborn \date 08/08
   */
  void Zero();

  //@}

protected:
  /*!
   * \brief Internal reshape function
   */
  int DoReshape(const int NumRows, const int NumCols, const bool light);

  std::size_t allocatedSize_; /// Remembers the amount of bytes allocated when changing layout

}; // class SerialDenseMatrix
} // namespace LINALG










#endif  // #ifndef SERIALDENSEMATRIX_H
