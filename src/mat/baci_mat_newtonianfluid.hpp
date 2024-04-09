/*----------------------------------------------------------------------*/
/*! \file
\brief Newtonian fluid material

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_NEWTONIANFLUID_HPP
#define FOUR_C_MAT_NEWTONIANFLUID_HPP



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
    /// material parameters for Newtonian fluid
    ///
    /// This object exists only once for each read Newton fluid.
    class NewtonianFluid : public Parameter
    {
     public:
      /// standard constructor
      NewtonianFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// kinematic or dynamic viscosity
      const double viscosity_;
      /// density
      const double density_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class NewtonianFluid

  }  // namespace PAR

  class NewtonianFluidType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "NewtonianFluidType"; }

    static NewtonianFluidType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static NewtonianFluidType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Newtonian fluid material
  ///
  /// This object exists (several times) at every element
  class NewtonianFluid : public Material
  {
   public:
    /// construct empty material object
    NewtonianFluid();

    /// construct the material object given material parameters
    explicit NewtonianFluid(MAT::PAR::NewtonianFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return NewtonianFluidType::Instance().UniqueParObjectId();
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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_fluid; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new NewtonianFluid(*this));
    }

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return density
    double Density() const override { return params_->density_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::NewtonianFluid* params_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
