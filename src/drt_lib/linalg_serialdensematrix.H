/*!----------------------------------------------------------------------
\file linalg_serialdensematrix.H

\class LINALG::SerialDenseMatrix

\brief A class that wraps Epetra_SerialDenseMatrix with minor modifications
       in the constructor

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef SERIALDENSEMATRIX_H
#define SERIALDENSEMATRIX_H


#include "Epetra_SerialDenseMatrix.h"


/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class SerialDenseMatrix

\brief A class that wraps Epetra_SerialDenseMatrix with minor modifications
       in the constructor

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SerialDenseMatrix : public Epetra_SerialDenseMatrix
{
public:

  //! @name Constructors and destructors and related methods

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

  */
  SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA, int NumRows, int NumCols,
                           bool set_object_label=true);

  /*!
  \brief Standard Copy Constructor wraps
         Epetra_SerialDenseMatrix(const SerialDenseMatrix& Source);

  */
  SerialDenseMatrix(const SerialDenseMatrix& Source);
  
  /*!
  \brief Shape dimensions of this matrix but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseMatrix::Shape(...))

  */
  int LightShape(int NumRows, int NumCols);

  /*!
  \brief Reshape dimensions of this matrix but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseMatrix::Reshape(...))
         
  Existing values are kept as far as they fir into the new shape but
  potentially excess space in the matrix is NOT zeroed out

  */
  int LightReshape(int NumRows, int NumCols);

  /*!
  \brief Destructor

  */
  virtual ~SerialDenseMatrix();


  //@}

protected:


}; // class SerialDenseMatrix
} // namespace LINALG










#endif  // #ifndef SERIALDENSEMATRIX_H
#endif  // #ifdef CCADISCRET
