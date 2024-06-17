/*----------------------------------------------------------------------*/
/*! \file
\brief
Former file of Ursula Mayer

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_CARREAUYASUDA_HPP
#define FOUR_C_MAT_CARREAUYASUDA_HPP



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
    /// material parameters
    class CarreauYasuda : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      CarreauYasuda(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      const double nu_0_;    /* zero-shear viscosity */
      const double nu_inf_;  /* infinite-shear viscosity */
      const double lambda_;  /* characteristic time */
      const double a_param_; /* constant parameter */
      const double b_param_; /* constant parameter */
      const double density_; /* density */

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class CarreauYasuda

  }  // namespace PAR

  class CarreauYasudaType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "CarreauYasudaType"; }

    static CarreauYasudaType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static CarreauYasudaType instance_;
  };

  /// Nonlinear viscosity according to Carreau-Yasuda
  class CarreauYasuda : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    CarreauYasuda();

    /// construct the material object given material parameters
    explicit CarreauYasuda(Mat::PAR::CarreauYasuda* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */

    int UniqueParObjectId() const override
    {
      return CarreauYasudaType::Instance().UniqueParObjectId();
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
      return Core::Materials::m_carreauyasuda;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new CarreauYasuda(*this));
    }

    /// return material parameters for element calculation
    //@{

    /// return parameter for zero-shear viscosity
    double Nu0() const { return params_->nu_0_; }
    /// return parameter for infinite-shear viscosity
    double NuInf() const { return params_->nu_inf_; }
    /// parameter for characteristic time
    double Lambda() const { return params_->lambda_; }
    /// constant parameter
    double AParam() const { return params_->a_param_; }
    /// constant parameter
    double BParam() const { return params_->b_param_; }
    /// density
    double Density() const override { return params_->density_; }

    //@}

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::CarreauYasuda* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
