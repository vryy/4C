/*----------------------------------------------------------------------*/
/*!
\file hemoglobin_0d_O2_saturation.H

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef HEMOGLOBIN_0D_O2_SATURATION_H
#define HEMOGLOBIN_0D_O2_SATURATION_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for Hemoglobin 0D O2 saturation material
    ///
    class Hemoglobin_0d_O2_saturation : public Parameter
    {
    public:

      /// standard constructor
      Hemoglobin_0d_O2_saturation(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Hemoglobin_0d_O2_saturation() { }

      /// @name material parameters
      //@{
      /// how much of blood satisfies this rule
      const double per_volume_blood_;
      /// saturation volume of O2 in blood (usually 21.36ml of O2 per 100ml of
      /// blood)
      const double o2_sat_per_vol_blood_;
      /// p1/2 indicates the PO2 at which blood is 50% saturated
      const double p_half_;
      /// is the power of the sigmoidal like function
      const double power_;
      /// number of O2 molecules per volume of O2
      const double nO2_per_VO2_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Hemoglobin_0d_O2_saturation

  }  // namespace PAR

  class Hemoglobin_0d_O2_saturationType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "hemoglobin_0d_O2_saturationType"; }

    static Hemoglobin_0d_O2_saturationType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static Hemoglobin_0d_O2_saturationType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Hemoglobin 0D O2 saturation material
  ///
  /// This object exists (several times) at every element
  class Hemoglobin_0d_O2_saturation : public Material
  {
  public:

    /// construct empty material object
    Hemoglobin_0d_O2_saturation();

    /// construct the material object given material parameters
    explicit Hemoglobin_0d_O2_saturation(MAT::PAR::Hemoglobin_0d_O2_saturation* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return Hemoglobin_0d_O2_saturationType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_0d_o2_hemoglobin_saturation; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Hemoglobin_0d_O2_saturation( *this ) ); }

    ///how much of blood satisfies this rule
    double PerVolumeBlood() const {return params_->per_volume_blood_;}
#if 0
    ///The dominating value of E at zero O2 saturation
    double EoLower() const {return params_->eo_lower_;}

    ///Lower decay constant
    double TauLower() const {return params_->tau_lower_;}

    ///Lower VO2 translator
    double VO2oLower() const {return params_->vO2o_lower_;}

    ///Translating E at maximum O2 saturation
    double Eo() const {return params_->eo_;}

    ///The dominating value of E at maximum O2 saturation
    double EoUpper() const {return params_->eo_upper_;}

    ///Upper decay constant
    double TauUpper() const {return params_->tau_upper_;}

    ///Upper VO2 translator
    double VO2oUpper() const {return params_->vO2o_upper_;}
#else
    double O2SaturationPerVolumeBlood() const {return params_->o2_sat_per_vol_blood_;}
    double ReferencePressure () const {return params_->p_half_;}
    double SigmoidalPower() const {return params_->power_;}
#endif
    ///Number of O2 moles per VO2
    double NumO2PerVO2() const { return params_->nO2_per_VO2_;}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::Hemoglobin_0d_O2_saturation* params_;
  };

}

#endif

