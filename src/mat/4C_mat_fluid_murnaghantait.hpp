/*----------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid according to Murnaghan-Tait

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_FLUID_MURNAGHANTAIT_HPP
#define FOUR_C_MAT_FLUID_MURNAGHANTAIT_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for weakly compressible fluid according to Murnaghan-Tait
    ///
    /// This object exists only once for each read Murnaghan-Tait fluid.
    class MurnaghanTaitFluid : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      MurnaghanTaitFluid(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// viscosity
      const double viscosity_;
      /// reference density
      const double refdensity_;
      /// reference pressure
      const double refpressure_;
      /// reference reference bulk modulus
      const double refbulkmodulus_;
      /// material parameter according to Murnaghan-Tait
      const double matparameter_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class MurnaghanTaitFluid

  }  // namespace PAR

  class MurnaghanTaitFluidType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "MurnaghanTaitFluidType"; }

    static MurnaghanTaitFluidType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MurnaghanTaitFluidType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Murnaghan-Tait fluid material
  ///
  /// This object exists (several times) at every element
  class MurnaghanTaitFluid : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    MurnaghanTaitFluid();

    /// construct the material object given material parameters
    explicit MurnaghanTaitFluid(Mat::PAR::MurnaghanTaitFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return MurnaghanTaitFluidType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(Core::Communication::PackBuffer& data) const override;

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
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_fluid_murnaghantait;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new MurnaghanTaitFluid(*this));
    }

    /// compute density
    double ComputeDensity(const double press) const;

    /// return material parameters for element calculation
    //@{

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return reference density
    double RefDensity() const { return params_->refdensity_; }

    /// return reference pressure
    double RefPressure() const { return params_->refpressure_; }

    /// return reference bulk modulus
    double RefBulkModulus() const { return params_->refbulkmodulus_; }

    /// return material parameter according to Murnaghan-Tait
    double MatParameter() const { return params_->matparameter_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::MurnaghanTaitFluid* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
