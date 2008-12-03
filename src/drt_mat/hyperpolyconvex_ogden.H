/*!----------------------------------------------------------------------
\file hyperpolyconvex_ogden.H
\brief
This file contains the routines required for isotropic nearly
incompressible soft tissue with particular application to alveolar
parenchyma.

The input line should read

MAT 1 MAT_Struct_HyperPolyconvexOgden YOUNG 8600. NUE 0.499 C 1000. K1 1500. K2 8.5 DENS 0.001

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.delete/Members/wiechert
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifndef HYPERPOLYOGDEN_H_
#define HYPERPOLYOGDEN_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  class HyperPolyOgden : public Material
  {

  public:

    // empty constructor
    HyperPolyOgden();

    // constructor with given material parameters
    HyperPolyOgden(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_HyperPolyOgden; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_hyperpolyogden; }

    // THE material routine
    void Evaluate(LINALG::Matrix<6,1>* glstrain,
		  LINALG::Matrix<6,6>* cmat,
		  LINALG::Matrix<6,1>* stress);

  private:

    /// my material parameters
    MATERIAL* matdata_;
  };
}

#endif  // CCADISCRET
#endif  // HYPERPOLYOGDEN_H_
