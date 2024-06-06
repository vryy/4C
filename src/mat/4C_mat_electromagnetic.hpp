/*----------------------------------------------------------------------*/
/*! \file
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       example input line:
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732

\level 2


 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ELECTROMAGNETIC_HPP
#define FOUR_C_MAT_ELECTROMAGNETIC_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for electromagnetic wave propagation
    class ElectromagneticMat : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      ElectromagneticMat(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      enum Matparamnames
      {
        sigma_,
        epsilon_,
        mu_,
        first = sigma_,
        last = mu_
      };

    };  // class ElectromagneticMat

  }  // namespace PAR

  class ElectromagneticMatType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ElectromagneticMatType"; }

    static ElectromagneticMatType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ElectromagneticMatType instance_;

  };  // class ElectromagneticMatType

  /*----------------------------------------------------------------------*/
  /// Wrapper for Sound propagation material
  class ElectromagneticMat : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    ElectromagneticMat();

    /// construct the material object given material parameters
    explicit ElectromagneticMat(Mat::PAR::ElectromagneticMat* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return ElectromagneticMatType::Instance().UniqueParObjectId();
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

    //! @name Access methods

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_electromagneticmat;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ElectromagneticMat(*this));
    }

    /// conductivity
    double sigma(int eleid = -1) const { return params_->GetParameter(params_->sigma_, eleid); }

    /// permittivity coefficient
    double epsilon(int eleid = -1) const { return params_->GetParameter(params_->epsilon_, eleid); }

    /// permeability coefficient
    double mu(int eleid = -1) const { return params_->GetParameter(params_->mu_, eleid); }



    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    //@}


   private:
    /// my material parameters
    Mat::PAR::ElectromagneticMat* params_;
  };

}  // namespace Mat



FOUR_C_NAMESPACE_CLOSE

#endif
