/*----------------------------------------------------------------------*/
/*! \file
\brief scatra_mat_aniso.H

\level 3

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SCATRA_ANISO_HPP
#define FOUR_C_MAT_SCATRA_ANISO_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material.hpp"
#include "4C_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for anisotropic scalar transport material
    class ScatraMatAniso : public Parameter
    {
     public:
      /// standard constructor
      ScatraMatAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      enum Matparamnames
      {
        diff1,
        diff2,
        diff3,
        first = diff1,
        last = diff3
      };

    };  // class Scatra

  }  // namespace PAR

  class ScatraMatAnisoType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ScatraMatAnisoType"; }

    static ScatraMatAnisoType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraMatAnisoType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for anisotropic scalar transport material
  class ScatraMatAniso : public Material
  {
   public:
    /// construct empty material object
    ScatraMatAniso();

    /// construct the material object given material parameters
    explicit ScatraMatAniso(MAT::PAR::ScatraMatAniso* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ScatraMatAnisoType::Instance().UniqueParObjectId();
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
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_scatra_aniso;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraMatAniso(*this));
    }

    /// diffusivity
    CORE::LINALG::Matrix<3, 1> Diffusivity(int eleid = -1) const
    {
      CORE::LINALG::Matrix<3, 1> diff;
      diff(0, 0) = params_->GetParameter(params_->diff1, eleid);
      diff(1, 0) = params_->GetParameter(params_->diff2, eleid);
      diff(2, 0) = params_->GetParameter(params_->diff3, eleid);

      return diff;
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatAniso* params_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
