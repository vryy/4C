/*!----------------------------------------------------------------------
\file aaaneohooke.H
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke YOUNG 1.044E7 BETA 188.1E5 DENS 1.0

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/
#ifndef AAANEOHOOKE_H_
#define AAANEOHOOKE_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  class AAAneohooke : public Material
  {

  public:

    // empty constructor
    AAAneohooke();

    // constructor with given material parameters
    AAAneohooke(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_AAAneohooke; }

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

    double Density();

    // material type
    virtual MATERIAL_TYP MaterialType() const { return m_aaaneohooke; }

    // THE material routine
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e,
		  Epetra_SerialDenseMatrix* cmat_e,
		  Epetra_SerialDenseVector* stress_e);

    // THE material routine
    void Evaluate(const LINALG::FixedSizeSerialDenseMatrix<6,1>& glstrain,
		  LINALG::FixedSizeSerialDenseMatrix<6,6>& cmat,
		  LINALG::FixedSizeSerialDenseMatrix<6,1>& stress);

  private:

    /// my material parameters
    MATERIAL* matdata_;
  };
}

#endif  // CCADISCRET
#endif  // AAANEOHOOKE_H_
