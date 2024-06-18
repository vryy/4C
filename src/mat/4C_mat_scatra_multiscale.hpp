/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SCATRA_MULTISCALE_HPP
#define FOUR_C_MAT_SCATRA_MULTISCALE_HPP

#include "4C_config.hpp"

#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_micro_macro_coupling.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  class ScatraMultiScaleGP;

  namespace PAR
  {
    //! material parameters
    class ScatraMultiScale : public ScatraMat, public ScatraMicroMacroCoupling
    {
     public:
      //! constructor
      ScatraMultiScale(const Core::Mat::PAR::Parameter::Data& matdata);


      //! create material
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      //! return porosity
      double Porosity() const { return porosity_; };

      //! return tortuosity
      double Tortuosity() const { return tortuosity_; };

     private:
      //! @name material parameters
      //@{
      //! porosity
      const double porosity_;

      //! tortuosity
      const double tortuosity_;
      //@}
    };  // class Mat::PAR::ScatraMultiScale
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ScatraMultiScaleType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ScatraMultiScaleType"; };

    static ScatraMultiScaleType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMultiScale : public ScatraMat, public ScatraMicroMacroCoupling
  {
   public:
    //! construct empty material
    ScatraMultiScale();

    //! construct material with specific material parameters
    explicit ScatraMultiScale(Mat::PAR::ScatraMultiScale* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ScatraMultiScaleType::Instance().UniqueParObjectId();
    };

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
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

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void unpack(const std::vector<char>& data) override;
    //@}

    //! return material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_scatra_multiscale;
    };

    //! clone material
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraMultiScale(*this));
    };

    //! return porosity
    double Porosity() const { return params_->Porosity(); };

    //! return tortuosity
    double Tortuosity() const { return params_->Tortuosity(); };

   private:
    //! return material parameters
    const Mat::PAR::ScatraMicroMacroCoupling* params() const override { return params_; };

    //! material parameters
    Mat::PAR::ScatraMultiScale* params_;
  };
}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
