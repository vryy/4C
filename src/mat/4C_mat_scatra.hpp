/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material

\level 1

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SCATRA_HPP
#define FOUR_C_MAT_SCATRA_HPP



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
    /// parameters for scalar transport material
    class ScatraMat : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      ScatraMat(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      enum Matparamnames
      {
        diff,
        reac,
        densific,
        /*!< True if the scalar reacts to the external force.
         * Example: the scalar is magnetic and reacts to the magnetic field. Another scalar could be
         * not magnetic and not react to the magnetic field. */
        reacts_to_external_force,
        first = diff,
        last = reacts_to_external_force
      };

    };  // class Scatra

  }  // namespace PAR

  class ScatraMatType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ScatraMatType"; }

    static ScatraMatType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMat : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    ScatraMat();

    /// construct the material object given material parameters
    explicit ScatraMat(Mat::PAR::ScatraMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return ScatraMatType::Instance().UniqueParObjectId(); }

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
      return Core::Materials::m_scatra;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraMat(*this));
    }

    /// diffusivity
    double Diffusivity(int eleid = -1) const { return params_->GetParameter(params_->diff, eleid); }

    /// reaction coefficient
    double ReaCoeff(int eleid = -1) const { return params_->GetParameter(params_->reac, eleid); }

    /// densification coefficient
    double Densification(int eleid = -1) const
    {
      return params_->GetParameter(params_->densific, eleid);
    }

    /// reacts to external force
    [[nodiscard]] double reacts_to_external_force() const
    {
      return params_->GetParameter(params_->reacts_to_external_force, -1);
    }


    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::ScatraMat* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
