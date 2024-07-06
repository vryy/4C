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
    std::string name() const override { return "IonType"; }

    static IonType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

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
    int unique_par_object_id() const override { return IonType::instance().unique_par_object_id(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by unique_par_object_id() which will then
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
      unique_par_object_id().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    Core::Materials::MaterialType material_type() const override { return Core::Materials::m_ion; }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> clone() const override
    {
      return Teuchos::rcp(new Ion(*this));
    }

    /// valence (= charge number)
    double valence() const { return params_->valence_; }
    /// diffusivity coefficient
    double diffusivity() const { return params_->diffusivity_; }
    /// densification coefficient
    double densification() const { return params_->densification_; }
    /// valence (= charge number) of eliminated ion species
    double elim_valence() const { return params_->elimvalence_; }
    /// diffusivity coefficient of eliminated ion species
    double elim_diffusivity() const { return params_->elimdiffusivity_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::Ion* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
