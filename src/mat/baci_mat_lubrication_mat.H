/*--------------------------------------------------------------------------*/
/*! \file
\brief Material model for the lubrication film

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_MAT_LUBRICATION_MAT_H
#define BACI_MAT_LUBRICATION_MAT_H


#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    class LubricationLaw;

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class LubricationMat : public Parameter
    {
     public:
      /// standard constructor
      LubricationMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// density
      const double density_;

      /// lubrication law ID
      int lubricationlawID_;

      //@}

      // implementation of lubrication law
      LubricationLaw* lubricationlaw_;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class Lubrication

  }  // namespace PAR

  class LubricationMatType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "LubricationMatType"; }

    static LubricationMatType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static LubricationMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class LubricationMat : public Material
  {
   public:
    /// construct empty material object
    LubricationMat();

    /// construct the material object given material parameters
    explicit LubricationMat(MAT::PAR::LubricationMat* params);

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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_lubrication; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new LubricationMat(*this));
    }

    /// compute current viscosity
    double ComputeViscosity(const double press  ///< (i) lubrication pressure
    );

    //! evaluate constitutive relation for viscosity and compute derivatives
    double ComputeViscosityDeriv(const double press, const double visc);

    /// density
    double Density() const override { return params_->density_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::LubricationMat* params_;
  };

}  // namespace MAT


BACI_NAMESPACE_CLOSE

#endif  // MAT_LUBRICATION_MAT_H
