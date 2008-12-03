/*!----------------------------------------------------------------------
\file lung_penalty.H
\brief

<pre>
Maintainer: Lena Wiechert & Sophie Rausch
            {wiechert,rausch}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*----------------------------------------------------------------------*/
#ifndef LUNGPENALTY_H_
#define LUNGPENALTY_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for hyperelastic polyconvex anisotropic material
  class LungPenalty : public Material
  {
  public:

    /// construct empty material object
    LungPenalty();

    /// construct the material object given material parameters
    explicit LungPenalty(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_LungPenalty; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_lung_penalty; }

    // wrapper for less intensive error output information
    void Evaluate(const LINALG::Matrix<6,1>* glstrain,
                  LINALG::Matrix<6,6>* cmat,
                  LINALG::Matrix<6,1>* stress);


    double Density();


  private:

    /// my material parameters
    MATERIAL* matdata_;

  };
}

#endif /*CCADISCRET*/
#endif /*LUNGPENALTY_H_*/
