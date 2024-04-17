/*----------------------------------------------------------------------*/
/*! \file
\brief material for heat transport due to Fourier-type thermal conduction and the Soret effect


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SORET_HPP
#define FOUR_C_MAT_SORET_HPP

#include "baci_config.hpp"

#include "baci_mat_fourieriso.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    //! parameters for Soret material
    class Soret : public FourierIso
    {
     public:
      //! constructor
      Soret(Teuchos::RCP<MAT::PAR::Material> matdata);


      //! create instance of Soret material
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! return Soret coefficient
      double SoretCoefficient() const { return soretcoefficient_; };

     private:
      //! @name parameters for Soret material
      //@{
      //! Soret coefficient
      const double soretcoefficient_;
      //@}
    };  // class MAT::PAR::Soret
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class SoretType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "SoretType"; };

    static SoretType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static SoretType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for Soret material
  class Soret : public FourierIso
  {
   public:
    //! construct empty Soret material
    Soret();

    //! construct Soret material with specific material parameters
    explicit Soret(MAT::PAR::Soret* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return SoretType::Instance().UniqueParObjectId(); };

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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_soret; };

    //! clone Soret material
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Soret(*this)); };

    //! return Soret coefficient
    double SoretCoefficient() const { return params_->SoretCoefficient(); };

   private:
    //! return material parameters
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    //! material parameters
    MAT::PAR::Soret* params_;
  };
}  // namespace MAT
FOUR_C_NAMESPACE_CLOSE

#endif
