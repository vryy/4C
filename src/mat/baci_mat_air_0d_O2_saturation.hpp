/*----------------------------------------------------------------------*/
/*! \file

\brief Gives relevant quantities of O2 saturation of air, used for scatra in reduced dimensional
airway elements framework (transport in elements and between air and blood)


\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_AIR_0D_O2_SATURATION_HPP
#define FOUR_C_MAT_AIR_0D_O2_SATURATION_HPP



#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_par_parameter.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Air 0D O2 saturation material
    ///
    class Air_0d_O2_saturation : public Parameter
    {
     public:
      /// standard constructor
      Air_0d_O2_saturation(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// Atmospheric pressure
      const double atmospheric_p_;
      /// number of O2 moles per volume unit of air
      const double nO2_per_VO2_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class Air_0d_O2_saturation

  }  // namespace PAR

  class Air_0d_O2_saturationType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "air_0d_O2_saturationType"; }

    static Air_0d_O2_saturationType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Air_0d_O2_saturationType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Air 0D O2 saturation material
  ///
  /// This object exists (several times) at every element
  class Air_0d_O2_saturation : public Material
  {
   public:
    /// construct empty material object
    Air_0d_O2_saturation();

    /// construct the material object given material parameters
    explicit Air_0d_O2_saturation(MAT::PAR::Air_0d_O2_saturation* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return Air_0d_O2_saturationType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_0d_o2_air_saturation;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Air_0d_O2_saturation(*this));
    }

    /// Atmospheric pressure
    double AtmosphericPressure() const { return params_->atmospheric_p_; }

    /// Number of O2 moles per unit volume of air
    double NumO2PerVO2() const { return params_->nO2_per_VO2_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Air_0d_O2_saturation* params_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
