/*--------------------------------------------------------------------------*/
/*! \file
\brief Material model for the lubrication film

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_LUBRICATION_MAT_HPP
#define FOUR_C_MAT_LUBRICATION_MAT_HPP


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
    class LubricationLaw;

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class LubricationMat : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      LubricationMat(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// density
      const double density_;

      /// lubrication law ID
      int lubricationlawID_;

      //@}

      // implementation of lubrication law
      LubricationLaw* lubricationlaw_;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class Lubrication

  }  // namespace PAR

  class LubricationMatType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "LubricationMatType"; }

    static LubricationMatType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LubricationMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class LubricationMat : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    LubricationMat();

    /// construct the material object given material parameters
    explicit LubricationMat(Mat::PAR::LubricationMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return LubricationMatType::Instance().UniqueParObjectId();
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
      return Core::Materials::m_lubrication;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new LubricationMat(*this));
    }

    /// compute current viscosity
    double ComputeViscosity(const double press  ///< (i) lubrication pressure
    );

    //! evaluate constitutive relation for viscosity and compute derivatives
    double compute_viscosity_deriv(const double press, const double visc);

    /// density
    double Density() const override { return params_->density_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::LubricationMat* params_;
  };

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif
