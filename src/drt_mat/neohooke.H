/*----------------------------------------------------------------------*/
/*!
\file neohooke.H
\brief contains the functions to establish local material law
       stress-strain law for isotropic material for a 3D hex element
       following compressible Neo Hookean material law
       see Holzapfel, Nonlinear Solid Mechanics, pp. 243
       example input line:
       MAT 1 MAT_Struct_NeoHooke  YOUNG 100.0 NUE 0.3 DENS 1.0
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef NEOHOOKE_H_
#define NEOHOOKE_H_

#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for Neo-Hookean Material
  class NeoHooke : public Material
  {
  public:

    /// construct empty material object
    NeoHooke();

    /// construct the material object given material parameters
    explicit NeoHooke(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_NeoHooke; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual MATERIAL_TYP MaterialType() const { return m_neohooke; }
    
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e,
                  Epetra_SerialDenseMatrix* cmat_e,
                  Epetra_SerialDenseVector* stress_e);
    
    void Evaluate(const LINALG::FixedSizeSerialDenseMatrix<6,1>& glstrain,
                  LINALG::FixedSizeSerialDenseMatrix<6,6>& cmat,
                  LINALG::FixedSizeSerialDenseMatrix<6,1>& stress);

    inline double Density() { return matdata_->m.neohooke->density; }


  private:

    /// my material parameters
    MATERIAL* matdata_;

    void InverseTensor(
                  const Epetra_SerialDenseMatrix& M,
                  Epetra_SerialDenseMatrix& Minv,
		  const double I3);

  
  };
}

#endif /*CCADISCRET*/
#endif /*NEOHOOKE_H_*/
