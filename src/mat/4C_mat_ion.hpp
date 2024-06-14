/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores parameters of ion species in electrolyte solution. Former file of Georg
Bauer

\level 2


*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ION_HPP
#define FOUR_C_MAT_ION_HPP

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
    /// material parameters for convection-diffusion
    class Ion : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      Ion(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// valence (= charge number)
      const double valence_;
      /// diffusivity coefficient
      const double diffusivity_;
      /// densification coefficient
      const double densification_;
      /// valence (= charge number) of eliminated ion species
      const double elimvalence_;
      /// diffusivity coefficient of eliminated ion species
      const double elimdiffusivity_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class Ion

  }  // namespace PAR

  class IonType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "IonType"; }

    static IonType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static IonType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Ion : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    Ion();

    /// construct the material object given material parameters
    explicit Ion(Mat::PAR::Ion* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return IonType::Instance().UniqueParObjectId(); }

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
    Core::Materials::MaterialType MaterialType() const override { return Core::Materials::m_ion; }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new Ion(*this));
    }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }
    /// diffusivity coefficient
    double Diffusivity() const { return params_->diffusivity_; }
    /// densification coefficient
    double Densification() const { return params_->densification_; }
    /// valence (= charge number) of eliminated ion species
    double ElimValence() const { return params_->elimvalence_; }
    /// diffusivity coefficient of eliminated ion species
    double ElimDiffusivity() const { return params_->elimdiffusivity_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::Ion* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
