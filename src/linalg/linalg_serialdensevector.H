/*----------------------------------------------------------------------*/
/*!
\file linalg_serialdensevector.H

\brief Specification of wrapper to a serial vector

<pre>
\level 0
\maintainer Martin Kronbichler
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

  Epetra_SerialDenseVector always allocates new memory upon call of Size() or
  Resize(). This class re-uses previously allocated memory as much as possible,
  similar to the implementation in std::vector. If the vector size is reduced
  or if the allocated memory is sufficient for the new vector size, the old
  data array is re-used. Re-allocation only happens in case more memory is
  necessary. This way, some memory might be wasted but it results in code that
  is generally considerably faster and avoids memory fragmentation. To force a
  new allocation with the correct amount of memory, call Size(0) first.
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
    SerialDenseVector(int Length, bool init = false);

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
    \brief Set the length of this vector and initialize to zero

    If the given size is zero, the memory will be freed. Otherwise,
    the algorithm keeps the old memory if it is sufficient and only sets the
    size bound differently.
    */
    int Size(int Length);

    /*!
    \brief Set length of this vector but do NOT initialize to zero
           (as opposed to Epetra_SerialDenseVector::Shape(...))

    */
    int LightSize(int Length);

    /*!
    \brief Resize length of this vector, keep old values, set excess memory to zero

    If the given size is zero, the memory will be freed. Otherwise,
    the algorithm keeps the old memory if it is sufficient and only sets the
    size bound differently.
    */
    int Resize(int Length);

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
    void Update(const double& ScalarB,     /*!< scale input vector with */
        const Epetra_SerialDenseVector& B, /*!< input vector */
        const double& ScalarThis           /*!< scale this vector with */
    );

    /*!
     * \brief Set vector components to zero,
     *        this = 0.0
     * \author a.ger \date 11/08
     */
    void Zero();

    //@}

   protected:
    int allocatedSize_;  /// Remembers the amount of bytes allocated when changing vector length

  };  // class SerialDenseVector
}  // namespace LINALG


// << operator
// ostream& operator << (ostream& os, const LINALG::SerialDenseVector& vector);



#endif  // #ifndef SERIALDENSEVECTOR_H
