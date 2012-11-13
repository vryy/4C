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

#include "material.H"
#include "matpar_parameter.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "../drt_lib/drt_parobjectfactory.H"


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

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class MooneyRivlin

  }  // namespace PAR

class MooneyRivlinType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "MooneyRivlinType"; }

  static MooneyRivlinType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static MooneyRivlinType instance_;
};

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
    virtual int UniqueParObjectId() const { return MooneyRivlinType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_mooneyrivlin; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new MooneyRivlin( *this ) ); }

    /*!
    \brief Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.

    The underlying strain-energy function is (expressed in invariants I1 and I2):

    W = c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + kappa * (J - 1)^2 +
        lambda / 4 * (J^2 - 1 - 2 ln(J))

    For references see Holzapfel p. 245, (Klinkel et al 2007).

    Parameters are c_1, c_2 and kappa as penalty factor to effectuate incompressibility
    (shear modulus mu = (c1 + c2) / 2)

    Note that this special form is stress-free in reference configuration if lambda == 0.0
    (See Holzapfel for discussion with d = 2 c_1 + 4 c_2)
    */
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress);

    /// return density
    virtual double Density() const { return params_->density_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::MooneyRivlin* params_;

  };
}

#endif /*MOONEYRIVLIN_H_*/
