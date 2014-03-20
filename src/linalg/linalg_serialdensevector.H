/*----------------------------------------------------------------------*/
/*!
\file linalg_serialdensevector.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SERIALDENSEVECTOR_H
#define SERIALDENSEVECTOR_H


#include "Epetra_SerialDenseVector.h"


namespace LINALG
{

/*!
\brief A class that wraps Epetra_SerialDenseVector with minor modifications
       in the constructor

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SerialDenseVector : public Epetra_SerialDenseVector
{
public:

  //! @name Constructors and destructors and related methods
  //@{

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
  SerialDenseVector(int Length, bool init=false);

  /*!
  \brief Standard Constructor wrapping
         Epetra_SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);


  */
  SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);

  /*!
  \brief Standard Constructor allowing view on Epetra_SerialDenseVector and wrapping
         SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);

  */
  SerialDenseVector(Epetra_SerialDenseVector& Source, Epetra_DataAccess CV);

  /*!
  \brief Standard Copy Constructor wraps
         Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source);

  */
  SerialDenseVector(const SerialDenseVector& Source);

  /*!
  \brief Epetra_SerialDenseVector Copy Constructor wraps
         Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source);

  */
  SerialDenseVector(const Epetra_SerialDenseVector& Source);

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

  //! @name Mathematical methods
  //@{

  /*!
  \brief Update vector components with scaled values of B,
         this = ScalarThis * this + ScalarB * B
  \author bborn \date 08/08
  */
  void Update(
    const double& ScalarB,  /*!< scale input vector with */
    const Epetra_SerialDenseVector& B,  /*!< input vector */
    const double& ScalarThis  /*!< scale this vector with */
  );

  /*!
   * \brief Set vector components to zero,
   *        this = 0.0
   * \author a.ger \date 11/08
   */
  void Zero();

  //@}

protected:


}; // class SerialDenseVector
} // namespace LINALG


// << operator
// ostream& operator << (ostream& os, const LINALG::SerialDenseVector& vector);









#endif  // #ifndef SERIALDENSEVECTOR_H
