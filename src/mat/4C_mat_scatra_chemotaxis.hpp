/*----------------------------------------------------------------------*/
/*! \file
 \brief
This file contains the base material for chemotactic scalars.

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_SCATRA_CHEMOTAXIS_HPP
#define FOUR_C_MAT_SCATRA_CHEMOTAXIS_HPP

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
    class ScatraChemotaxisMat : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      ScatraChemotaxisMat(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// number of scalars in this pair
      const int numscal_;

      /// the list of material IDs
      const std::vector<int> pair_;

      /// reaction coefficient
      const double chemocoeff_;

    };  // class Scatra

  }  // namespace PAR

  class ScatraChemotaxisMatType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ScatraChemotaxisMatType"; }

    static ScatraChemotaxisMatType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScatraChemotaxisMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraChemotaxisMat : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    ScatraChemotaxisMat();

    /// construct the material object given material parameters
    explicit ScatraChemotaxisMat(Mat::PAR::ScatraChemotaxisMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ScatraChemotaxisMatType::Instance().UniqueParObjectId();
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
      return Core::Materials::m_scatra_chemotaxis;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ScatraChemotaxisMat(*this));
    }

    /// return number of scalars for this reaction
    int NumScal() const { return params_->numscal_; }

    /// return chemotaxis coefficient
    double ChemoCoeff() const { return params_->chemocoeff_; }

    /// return pairing
    const std::vector<int>* Pair() const { return &params_->pair_; }

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::ScatraChemotaxisMat* params_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
