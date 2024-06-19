/*----------------------------------------------------------------------*/
/*! \file
\brief temperature-dependent gas according to Sutherland law

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SUTHERLAND_HPP
#define FOUR_C_MAT_SUTHERLAND_HPP



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
    /// parameters for material with temperature dependence according to Sutherland law
    class Sutherland : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      Sutherland(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// reference dynamic viscosity (kg/(m*s))
      const double refvisc_;
      /// reference temperature (K)
      const double reftemp_;
      /// Sutherland temperature (K)
      const double suthtemp_;
      /// specific heat capacity
      const double shc_;
      /// Prandtl number
      const double pranum_;
      /// (initial) thermodynamic pressure
      const double thermpress_;
      /// specific gas constant R
      const double gasconst_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class Sutherland

  }  // namespace PAR

  class SutherlandType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "SutherlandType"; }

    static SutherlandType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static SutherlandType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// material with temperature dependence according to Sutherland law
  class Sutherland : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    Sutherland();

    /// construct the material object given material parameters
    explicit Sutherland(Mat::PAR::Sutherland* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */

    int UniqueParObjectId() const override
    {
      return SutherlandType::Instance().UniqueParObjectId();
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
      return Core::Materials::m_sutherland;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new Sutherland(*this));
    }

    /// compute viscosity
    double ComputeViscosity(const double temp) const;

    /// compute diffusivity
    double ComputeDiffusivity(const double temp) const;

    /// compute density
    double ComputeDensity(const double temp, const double thermpress) const;

    /// return material parameters for element calculation
    //@{

    /// reference dynamic viscosity (kg/(m*s))
    double RefVisc() const { return params_->refvisc_; }
    /// reference temperature (K)
    double RefTemp() const { return params_->reftemp_; }
    /// Sutherland temperature (K)
    double SuthTemp() const { return params_->suthtemp_; }
    /// specific heat capacity
    double Shc() const { return params_->shc_; }
    /// Prandtl number
    double PraNum() const { return params_->pranum_; }
    /// (initial) thermodynamic pressure
    double ThermPress() const { return params_->thermpress_; }
    /// specific gas constant R
    double GasConst() const { return params_->gasconst_; }

    //@}

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::Sutherland* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
