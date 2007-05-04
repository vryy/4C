/*!----------------------------------------------------------------------
\file linalg_serialdensevector.H

\class LINALG::SerialDenseVector

\brief A class that wraps Epetra_SerialDenseVector with minor modifications
       in the constructor

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SERIALDENSEVECTOR_H
#define SERIALDENSEVECTOR_H


#include "Epetra_SerialDenseVector.h"


/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class SerialDenseVector

\brief A class that wraps Epetra_SerialDenseVector with minor modifications
       in the constructor

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SerialDenseVector : public Epetra_SerialDenseVector
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseVector();
         
  */
  SerialDenseVector();

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseVector(int Length);
         
  Opposed to the base class constructor, this one does NOT initialize the matrix to zero.
  default behavior of set_object_label is also reversed compared to base class 
  Epetra_SerialDenseVector
  
  \param Length           (in): supposed length of vector

  */
  SerialDenseVector(int Length);

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);


  */
  SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);

  /*!
  \brief Standard Copy Constructor wraps
         Epetra_SerialDenseVector(const SerialDenseVector& Source);

  */
  SerialDenseVector(const SerialDenseVector& Source);
  
  /*!
  \brief Set length of this vector but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseVector::Shape(...))

  */
  int LightSize(int Length);

  /*!
  \brief Resize length of this vector but do NOT initialize to zero
         (as opposed to Epetra_SerialDenseVector::Reshape(...))
         
  Existing values are kept as far as they fir into the new shape but
  potentially excess space in the matrix is NOT zeroed out

  */
  int LightResize(int Length);

  /*!
  \brief Destructor

  */
  virtual ~SerialDenseVector();


  //@}

protected:


}; // class SerialDenseVector
} // namespace LINALG


// << operator
ostream& operator << (ostream& os, const LINALG::SerialDenseVector& vector);









#endif  // #ifndef SERIALDENSEVECTOR_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
