/*----------------------------------------------------------------------*/
/*! \file
\brief Linear law (pressure-dependent) for the density and the viscosity

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_FLUID_LINEAR_DENSITY_VISCOSITY_HPP
#define FOUR_C_MAT_FLUID_LINEAR_DENSITY_VISCOSITY_HPP



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
    /// material parameters for fluid with linear law (pressure-dependent) for
    /// the density and the viscosity
    ///
    /// This object exists only once for each read fluid.
    class LinearDensityViscosity : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      LinearDensityViscosity(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// reference density
      const double refdensity_;
      /// reference viscosity
      const double refviscosity_;
      /// reference pressure
      const double refpressure_;
      /// density-pressure coefficient
      const double coeffdensity_;
      /// viscosity-pressure coefficient
      const double coeffviscosity_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class LinearDensityViscosity

  }  // namespace PAR

  class LinearDensityViscosityType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "LinearDensityViscosityType"; }

    static LinearDensityViscosityType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LinearDensityViscosityType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for fluid with linear law (pressure-dependent) for
  /// the density and the viscosity
  ///
  /// This object exists (several times) at every element
  class LinearDensityViscosity : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    LinearDensityViscosity();

    /// construct the material object given material parameters
    explicit LinearDensityViscosity(Mat::PAR::LinearDensityViscosity* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return LinearDensityViscosityType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void pack(Core::Communication::PackBuffer& data) const override;

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
    void unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_fluid_linear_density_viscosity;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new LinearDensityViscosity(*this));
    }

    /// compute density
    double ComputeDensity(const double press) const;

    /// compute viscosity
    double ComputeViscosity(const double press) const;

    /// return material parameters for element calculation
    //@{

    /// return reference density
    double RefDensity() const { return params_->refdensity_; }

    /// return reference viscosity
    double RefViscosity() const { return params_->refviscosity_; }

    /// return reference pressure
    double RefPressure() const { return params_->refpressure_; }

    /// return density-pressure coefficient
    double CoeffDensity() const { return params_->coeffdensity_; }

    /// return viscosity-pressure coefficient
    double CoeffViscosity() const { return params_->coeffviscosity_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::LinearDensityViscosity* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
