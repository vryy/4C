/*----------------------------------------------------------------------*/
/*! \file
\brief Material law for elastic spring (either translational or rotational spring)

\level 3

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SPRING_HPP
#define FOUR_C_MAT_SPRING_HPP



#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_par_parameter.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for spring
    class Spring : public Parameter
    {
     public:
      /// standard constructor
      Spring(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// stiffness (translational or rotational)
      const double stiffness_;

      /// density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class Spring

  }  // namespace PAR

  class SpringType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "SpringType"; }

    static SpringType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static SpringType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for spring material
  class Spring : public Material
  {
   public:
    /// construct empty material object
    Spring();

    /// construct the material object given material parameters
    explicit Spring(MAT::PAR::Spring* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return SpringType::Instance().UniqueParObjectId(); }

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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_spring; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Spring(*this)); }

    /// stiffness (translational or rotational)
    double Stiffness() const { return params_->stiffness_; }

    /// density
    double Density() const override { return params_->density_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Spring* params_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
