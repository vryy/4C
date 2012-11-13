/*!----------------------------------------------------------------------
\file anisotropic_balzani.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef ANISOTROPIC_BALZANI_H_
#define ANISOTROPIC_BALZANI_H_

#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RCP.hpp>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for hyperelastic polyconvex anisotropic material
    class AnisotropicBalzani : public Parameter
    {
    public:

      /// standard constructor
      AnisotropicBalzani(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~AnisotropicBalzani() { ; }

      /// @name material parameters
      //@{

      const double c1_;
      const double eps1_;
      const double eps2_;
      const double alpha1_;
      const double alpha2_;
      const double dens_;
      const int aloc_;
      double a1_[3];
      const double alpha1_2_;
      const double alpha2_2_;
      double a2_[3];

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class AnisotropicBalzani

  }  // namespace PAR

class AnisotropicBalzaniType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "AnisotropicBalzaniType"; }

  static AnisotropicBalzaniType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static AnisotropicBalzaniType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for hyperelastic polyconvex anisotropic material
  class AnisotropicBalzani : public Material
  {
  public:

    /// construct empty material object
    AnisotropicBalzani();

    /// construct the material object given material parameters
    explicit AnisotropicBalzani(MAT::PAR::AnisotropicBalzani* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return AnisotropicBalzaniType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_anisotropic_balzani; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new AnisotropicBalzani( *this ) ); }

    /// evaluate stress and elasticity matrix
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            const int gp, const int ele_ID, const double time,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress);

    /// return density
    virtual double Density() const { return params_->dens_; }

    /*!
      \brief setup material at element input with potentially local fiber cosy

      This anisotropic material consists of fiber angles which can be defined
      at element level via local coordinate systems. This method controls setup
      of corresponding variables and reads e.g. fiber directions into its history
      variables.
    */
    void Setup(DRT::INPUT::LineDefinition* linedef);

    /// return fiber directions
    std::vector<double> ReturnGlobalFiberDirection() {
      std::vector<double> aglob(3);
      aglob[0] = params_->a1_[0];
      aglob[1] = params_->a1_[1];
      aglob[2] = params_->a1_[2];
      return aglob;
    }

    /// check whether fiber directions are globally defined
    bool GlobalFiberDirection() const {
      if (params_->aloc_ == 1) return false;
      else return true;
    };

    /// return a1
    std::vector<double> Geta1() const {return a1_;}

    /// return a2
    std::vector<double> Geta2() const {return a2_;}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::AnisotropicBalzani* params_;

    /// first fiber vector
    std::vector<double> a1_;
    /// second fiber vector
    std::vector<double> a2_;
  };
}

#endif /*ANISOTROPIC_BALZANI_H_*/
