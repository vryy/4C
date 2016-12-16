/*----------------------------------------------------------------------*/
/*!
 \file scatra_growth_scd.H

 \brief  scatra material used for scalar dependent growth modeling

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364

\level 3
 *----------------------------------------------------------------------*/


#ifndef SCATRA_GROWTH_SCD_H_
#define SCATRA_GROWTH_SCD_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraGrowthScd : public Parameter
    {
    public:

      /// standard constructor
      ScatraGrowthScd(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ScatraGrowthScd() { ; }

      /// @name material parameters
      //@{

      /// kinematic diffusivity
     const double diffusivity_;
     /// densitiy of structure material
     const double strdensity_;
     /// substrate uptake rate coefficient
     const double rearate_;
     /// substrate saturation coefficient
     const double satcoeff_;
      //@}

     /// source mass term
     const std::string* sourcemass_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ScatraGrowthScd

  } // namespace PAR

class ScatraGrowthScdType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ScatraGrowthScdType"; }

  static ScatraGrowthScdType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ScatraGrowthScdType instance_;
};

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraGrowthScd : public Material
  {
  public:

    /// construct empty material object
    ScatraGrowthScd();

    /// construct the material object given material parameters
    explicit ScatraGrowthScd(MAT::PAR::ScatraGrowthScd* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ScatraGrowthScdType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_growth_scd; }

    /// compute reaction coefficient
    double ComputeReactionCoeff(const double csnp, const double theta, const double dtheta, const double detFe) const;

    /// compute reaction coefficient derivative
    double ComputeReactionCoeffDeriv(const double csnp, const double theta, const double thetaold, const double dt) const;

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ScatraGrowthScd( *this ) ); }

    /// diffusivity
    double Diffusivity() const { return params_->diffusivity_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ScatraGrowthScd* params_;

  };

}


#endif /* SCATRA_GROWTH_SCD_H_ */
