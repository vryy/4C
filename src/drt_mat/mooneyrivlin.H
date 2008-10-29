/*!----------------------------------------------------------------------
\file mooneyrivlin.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef MOONEYRIVLIN_H_
#define MOONEYRIVLIN_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for Mooney-Rivlin material
  class MooneyRivlin : public Material
  {
  public:

    /// construct empty material object
    MooneyRivlin();

    /// construct the material object given material parameters
    explicit MooneyRivlin(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_MooneyRivlin; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_mooneyrivlin; }

    /*!
    \brief Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.

    The underlying strain-energy function is ('a' is 'alpha' and J = I3^(1/2)):

    W = sum_(p=1)^2 ( mu_p / a_p (l1^a_p + l2^a_p + l3^a_p) - mu_p ln(J) )
      + Lambda/4 * (J^2 - 1 - 2 ln(J) )

    which can be expressed in terms of invariants I1 and I2 as:

    W = mu_1/a_1 (I1 - 3)  +  mu_2/a_2 (I2 - 3)  -  mu_1 ln(J)  -  mu_2 ln(J)
      + Lambda/4 * (J^2 - 1 - 2 ln(J) )

    For references see Holzapfel p. 245, Simo&Miehe 1992, Klinkel et al 2007.

    Parameters are mu_1, mu_2, a_1, a_2 and Lambda = kappa - 2/3 mu as penalty factor
    to effectuate incompressibility (shear modulus mu = (mu_1 a_1 + mu_2 a_2)/2 )

    Mind that it is not stress-free in reference configuration!
    */
    void Evaluate(
            const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1>* glstrain,
            LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> * stress);

    /// return density
    double Density();

  private:

    /// my material parameters
    MATERIAL* matdata_;

  };
}

#endif /*CCADISCRET*/
#endif /*MOONEYRIVLIN_H_*/
