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

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class MooneyRivlin : public Parameter
    {
    public:

      /// standard constructor
      MooneyRivlin(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~MooneyRivlin() { ; }

      /// @name material parameters
      //@{
      const double c1_;
      const double c2_;
      const double kap_;
      const double lambda_;
      const double density_;
      //@}

    };  // class MooneyRivlin

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /// Wrapper for Mooney-Rivlin material
  class MooneyRivlin : public Material
  {
  public:

    /// construct empty material object
    MooneyRivlin();

    /// construct the material object given material parameters
    explicit MooneyRivlin(MAT::PAR::MooneyRivlin* params);

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_mooneyrivlin; }

    /*!
    \brief Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.

    The underlying strain-energy function is (expressed in invariants I1 and I2):
    
    W = c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + kappa * (J - 1)^2 +
        lambda / 4 * (J^2 - 1 - 2 ln(J))

    For references see Holzapfel p. 245, (Klinkel et al 2007).

    Parameters are c_1, c_2 and kappa as penalty factor to effectuate incompressibility 
    (shear modulus mu = (c1 + c2) / 2)

    Note that this special form is stress-free in reference configuration iff lambda == 0.0 
    (See Holzapfel for discussion with d = 2 c_1 + 4 c_2)
    */
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress);

    /// return density
    double Density() const { return params_->density_; }

  private:

    /// my material parameters
    MAT::PAR::MooneyRivlin* params_;

  };
}

#endif /*CCADISCRET*/
#endif /*MOONEYRIVLIN_H_*/
