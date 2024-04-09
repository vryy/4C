/*----------------------------------------------------------------------*/
/*! \file
\brief non-Newtonian fluid of Herschel-Bulkley type

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_HERSCHELBULKLEY_HPP
#define FOUR_C_MAT_HERSCHELBULKLEY_HPP



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
    /// material parameters
    class HerschelBulkley : public Parameter
    {
     public:
      /// standard constructor
      HerschelBulkley(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      const double tau0_;           /* yield stress */
      const double kfac_;           /* constant factor */
      const double nexp_;           /* exponent */
      const double mexp_;           /* exponent */
      const double lolimshearrate_; /* lower limit of shear rate */
      const double uplimshearrate_; /* upper limit of shear rate */
      const double density_;        /* density */

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class HerschelBulkley

  }  // namespace PAR

  class HerschelBulkleyType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "HerschelBulkleyType"; }

    static HerschelBulkleyType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static HerschelBulkleyType instance_;
  };

  /// Nonlinear viscosity according to HerschelBulkley
  class HerschelBulkley : public Material
  {
   public:
    /// construct empty material object
    HerschelBulkley();

    /// construct the material object given material parameters
    explicit HerschelBulkley(MAT::PAR::HerschelBulkley* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */

    int UniqueParObjectId() const override
    {
      return HerschelBulkleyType::Instance().UniqueParObjectId();
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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_herschelbulkley; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new HerschelBulkley(*this));
    }

    /// return material parameters for element calculation
    //@{

    /// yield stress
    double Tau0() const { return params_->tau0_; }
    /// constant factor
    double KFac() const { return params_->kfac_; }
    /// exponent
    double NExp() const { return params_->nexp_; }
    /// exponent
    double MExp() const { return params_->mexp_; }
    /// lower limit of shear rate
    double LoLimShearRate() const { return params_->lolimshearrate_; }
    /// upper limit of shear rate
    double UpLimShearRate() const { return params_->uplimshearrate_; }
    /// density
    double Density() const override { return params_->density_; }

    //@}

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::HerschelBulkley* params_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
