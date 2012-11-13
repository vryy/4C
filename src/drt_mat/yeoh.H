/*!----------------------------------------------------------------------
\file yeoh.H
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/
#ifndef YEOH_H_
#define YEOH_H_

#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class Yeoh : public Parameter
    {
    public:

      /// standard constructor
      Yeoh(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Yeoh() { ; }

      /// @name material parameters
      //@{
      const double c1_;
      const double c2_;
      const double c3_;
      const double kap_;
      const double density_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Yeoh

  }  // namespace PAR

class YeohType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "YeohType"; }

  static YeohType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static YeohType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for Yeoh material
  class Yeoh : public Material
  {
  public:

    /// construct empty material object
    Yeoh();

    /// construct the material object given material parameters
    explicit Yeoh(MAT::PAR::Yeoh* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return YeohType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_yeoh; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Yeoh( *this ) ); }

    /*!
    \brief Yeoh type nearly incompressible, hyperelastic 3D material law.

    The underlying strain-energy function in decoupled form is given by:

    W = c_1 (I1_bar - 3)  +  c_2 (I1_bar - 3)^2  +  c_3 (I1_bar-3)^3 + kappa * (J - 1)^2,
    where I1_bar is the first invariant of the modified Cauchy-Green tensor C_bar = J^{-2/3} C.


    For references see Holzapfel p. 248f, 255, 263

    */
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress);

    /// return density
    virtual double Density() const { return params_->density_; }

    /// Return shear modulus
    double ShearMod() const {return params_->c1_/2.0;};

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::Yeoh* params_;

  };
}

#endif /*YEOH_H_*/
