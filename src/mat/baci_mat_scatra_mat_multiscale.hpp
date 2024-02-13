/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_SCATRA_MAT_MULTISCALE_HPP
#define BACI_MAT_SCATRA_MAT_MULTISCALE_HPP

#include "baci_config.hpp"

#include "baci_mat_scatra_mat.hpp"
#include "baci_mat_scatra_multiscale.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  // forward declaration
  class ScatraMultiScaleGP;

  namespace PAR
  {
    //! material parameters
    class ScatraMatMultiScale : public ScatraMat, public ScatraMultiScale
    {
     public:
      //! constructor
      ScatraMatMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata);


      //! create material
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

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
    };  // class MAT::PAR::ScatraMatMultiScale
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ScatraMatMultiScaleType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ScatraMatMultiScaleType"; };

    static ScatraMatMultiScaleType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraMatMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMatMultiScale : public ScatraMat, public ScatraMultiScale
  {
   public:
    //! construct empty material
    ScatraMatMultiScale();

    //! construct material with specific material parameters
    explicit ScatraMatMultiScale(MAT::PAR::ScatraMatMultiScale* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ScatraMatMultiScaleType::Instance().UniqueParObjectId();
    };

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
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

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void Unpack(const std::vector<char>& data) override;
    //@}

    //! return material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_scatra_multiscale;
    };

    //! clone material
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraMatMultiScale(*this));
    };

    //! return porosity
    double Porosity() const { return params_->Porosity(); };

    //! return tortuosity
    double Tortuosity() const { return params_->Tortuosity(); };

   private:
    //! return material parameters
    const MAT::PAR::ScatraMultiScale* Params() const override { return params_; };

    //! material parameters
    MAT::PAR::ScatraMatMultiScale* params_;
  };
}  // namespace MAT
BACI_NAMESPACE_CLOSE

#endif
