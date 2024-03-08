/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of electrochemistry problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_NEWMAN_MULTISCALE_HPP
#define BACI_MAT_NEWMAN_MULTISCALE_HPP

#include "baci_config.hpp"

#include "baci_mat_newman.hpp"
#include "baci_mat_scatra_micro_macro_coupling.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    //! material parameters
    class NewmanMultiScale : public Newman, public ScatraMicroMacroCoupling
    {
     public:
      //! constructor
      NewmanMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata);


      //! create instance of Newman multi-scale material
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! return electronic conductivity
      double Sigma() const { return sigma_; };

     private:
      //! @name parameters for Newman multi-scale material
      //@{
      //! electronic conductivity
      const double sigma_;
      //@}
    };  // class MAT::PAR::NewmanMultiScale
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class NewmanMultiScaleType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "NewmanMultiScaleType"; };

    static NewmanMultiScaleType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static NewmanMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for Newman multi-scale material
  class NewmanMultiScale : public Newman, public ScatraMicroMacroCoupling
  {
   public:
    //! construct empty Newman multi-scale material
    NewmanMultiScale();

    //! construct Newman multi-scale material with specific material parameters
    explicit NewmanMultiScale(MAT::PAR::NewmanMultiScale* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return NewmanMultiScaleType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_newman_multiscale;
    };

    //! clone Newman multi-scale material
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new NewmanMultiScale(*this));
    };

    //! return electronic conductivity
    double Sigma() const { return params_->Sigma(); };

   private:
    //! return material parameters
    const MAT::PAR::ScatraMicroMacroCoupling* Params() const override { return params_; };

    //! material parameters
    MAT::PAR::NewmanMultiScale* params_;
  };  // wrapper for Newman multi-scale material
}  // namespace MAT
BACI_NAMESPACE_CLOSE

#endif
