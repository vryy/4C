/*----------------------------------------------------------------------*/
/*! \file
 \brief scatra material for transport within porous model with special implementations
        for ECM model


\level 3
 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_SCATRA_MAT_PORO_ECM_HPP
#define BACI_MAT_SCATRA_MAT_PORO_ECM_HPP

#include "baci_config.hpp"

#include "baci_mat_scatra_reaction_mat.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraMatPoroECM : public ScatraReactionMat
    {
     public:
      /// standard constructor
      ScatraMatPoroECM(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      double reacscale_;
    };
    // class Scatra

  }  // namespace PAR

  class ScatraMatPoroECMType : public ScatraReactionMatType
  {
   public:
    std::string Name() const override { return "ScatraMatPoroECMType"; }

    static ScatraMatPoroECMType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraMatPoroECMType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMatPoroECM : public ScatraReactionMat
  {
   public:
    /// construct empty material object
    ScatraMatPoroECM();

    /// construct the material object given material parameters
    explicit ScatraMatPoroECM(MAT::PAR::ScatraMatPoroECM* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return ScatraReactionMatType::Instance().UniqueParObjectId();
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
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_scatra_reaction_poroECM;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraMatPoroECM(*this));
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// return reaction coefficient
    virtual double ReacCoeff() const { return reaccoeff_; }

    /// compute reaction coefficient from structure chemical potential
    void ComputeReacCoeff(double chempot);

   private:
    /// my material parameters
    MAT::PAR::ScatraMatPoroECM* params_;

    /// reaction coefficient
    double reaccoeff_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_SCATRA_MAT_PORO_ECM_H
