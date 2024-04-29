/*---------------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid

\level 1

*/
/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_FLUID_WEAKLY_COMPRESSIBLE_HPP
#define FOUR_C_MAT_FLUID_WEAKLY_COMPRESSIBLE_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material.hpp"
#include "4C_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for weakly compressible fluid
    ///
    /// This object exists only once for each read fluid.
    class WeaklyCompressibleFluid : public Parameter
    {
     public:
      /// standard constructor
      WeaklyCompressibleFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// viscosity
      const double viscosity_;
      /// reference density
      const double refdensity_;
      /// reference pressure
      const double refpressure_;
      /// compressibility coefficient
      const double comprcoeff_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class WeaklyCompressibleFluid

  }  // namespace PAR

  class WeaklyCompressibleFluidType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "WeaklyCompressibleFluidType"; }

    static WeaklyCompressibleFluidType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static WeaklyCompressibleFluidType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for weakly compressible fluid
  ///
  /// This object exists (several times) at every element
  class WeaklyCompressibleFluid : public Material
  {
   public:
    /// construct empty material object
    WeaklyCompressibleFluid();

    /// construct the material object given material parameters
    explicit WeaklyCompressibleFluid(MAT::PAR::WeaklyCompressibleFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return WeaklyCompressibleFluidType::Instance().UniqueParObjectId();
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
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_fluid_weakly_compressible;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new WeaklyCompressibleFluid(*this));
    }

    /// compute density
    double ComputeDensity(const double press) const;

    /// compute pressure
    double ComputePressure(const double dens) const;

    /// return material parameters for element calculation
    //@{

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return reference density
    double RefDensity() const { return params_->refdensity_; }

    /// return reference pressure
    double RefPressure() const { return params_->refpressure_; }

    /// return compressibility coefficient
    double ComprCoeff() const { return params_->comprcoeff_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::WeaklyCompressibleFluid* params_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
