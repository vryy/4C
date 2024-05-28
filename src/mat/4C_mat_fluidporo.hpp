/*----------------------------------------------------------------------*/
/*! \file
\brief  fluid material for poroelasticity problems


\level 2
 *-----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_FLUIDPORO_HPP
#define FOUR_C_MAT_FLUIDPORO_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace FLUIDPORO
  {
    class PoroAnisotropyStrategyBase;
  }
  namespace PAR
  {
    enum PoroFlowType
    {
      undefined,
      darcy,
      darcy_brinkman
    };

    enum PoroFlowPermeabilityFunction
    {
      pf_undefined,
      constant,
      kozeny_carman,
      const_material_transverse,
      const_material_orthotropic,
      const_material_nodal_orthotropic
    };

    /*----------------------------------------------------------------------*/
    //! material parameters for fluid poro
    //! This object exists only once for each read Newton fluid.
    class FluidPoro : public CORE::MAT::PAR::Parameter
    {
     public:
      //! standard constructor
      FluidPoro(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      //! set initial porosity from structural material and calculate
      //! permeability_correction_factor_
      void SetInitialPorosity(double initial_porosity);

      //! @name material parameters
      //!@{

      //! kinematic or dynamic viscosity
      const double viscosity_;
      //! density
      const double density_;
      //! permeability
      const double permeability_;
      //! axial permeability for material transverse isotropy
      const double axial_permeability_;
      //! vector of orthotropic permeabilities
      std::vector<double> orthotropic_permeabilities_;
      //! flow type: Darcy or Darcy-Brinkman
      PoroFlowType type_;
      //! flag indicating varying permeability
      const bool varying_permeability_;
      //! permeability function type
      PoroFlowPermeabilityFunction permeability_func_;
      //! a correction factor to ensure the permeability set in the input file
      double permeability_correction_factor_;
      //! initial porosity
      double initial_porosity_;

      //! create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;
    };

  }  // namespace PAR

  class FluidPoroType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "FluidPoroType"; }

    static FluidPoroType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FluidPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  //! Wrapper for poro fluid flow material
  //! This object exists (several times) at every element
  class FluidPoro : public CORE::MAT::Material
  {
   public:
    //! construct empty material object
    FluidPoro();

    //! construct the material object given material parameters
    explicit FluidPoro(MAT::PAR::FluidPoro* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override { return FluidPoroType::Instance().UniqueParObjectId(); }

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

    //!@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_fluidporo;
    }

    //! return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new FluidPoro(*this));
    }

    //! compute reaction tensor - 2D
    void compute_reaction_tensor(CORE::LINALG::Matrix<2, 2>& reaction_tensor, const double& J,
        const double& porosity,
        const std::vector<std::vector<double>>& anisotropic_permeability_directions = {},
        const std::vector<double>& anisotropic_permeability_coeffs = {}) const;

    //! compute reaction tensor - 3D
    void compute_reaction_tensor(CORE::LINALG::Matrix<3, 3>& reaction_tensor, const double& J,
        const double& porosity,
        const std::vector<std::vector<double>>& anisotropic_permeability_directions = {},
        const std::vector<double>& anisotropic_permeability_coeffs = {}) const;

    //! compute reaction coefficient
    double compute_reaction_coeff() const;

    //! compute linearization of reaction tensor - 2D
    void compute_lin_mat_reaction_tensor(CORE::LINALG::Matrix<2, 2>& linreac_dphi,
        CORE::LINALG::Matrix<2, 2>& linreac_dJ, const double& J, const double& porosity) const;

    //! compute linearization of reaction tensor - 3D
    void compute_lin_mat_reaction_tensor(CORE::LINALG::Matrix<3, 3>& linreac_dphi,
        CORE::LINALG::Matrix<3, 3>& linreac_dJ, const double& J, const double& porosity) const;

    //! effective viscosity (zero for Darcy and greater than zero for Darcy-Brinkman)
    double EffectiveViscosity() const;

    //! return type
    PAR::PoroFlowType Type() const { return params_->type_; }

    //! return viscosity
    double Viscosity() const { return params_->viscosity_; }

    //! return density
    double Density() const override { return params_->density_; }

    //! return permeability function
    PAR::PoroFlowPermeabilityFunction permeability_function() const
    {
      return params_->permeability_func_;
    }

    //! Return quick accessible material parameter data
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

    //! flag indicating a varying permeability
    bool VaryingPermeability() const { return params_->varying_permeability_; }

    //! flag indicating nodal orthotropy
    bool IsNodalOrthotropic() const
    {
      return params_->permeability_func_ == MAT::PAR::const_material_nodal_orthotropic;
    }

   private:
    //! my material parameters
    MAT::PAR::FluidPoro* params_;

    //! anisotropy strategy (isotropy, transverse isotropy, orthotropy, or nodal orthotropy)
    Teuchos::RCP<MAT::FLUIDPORO::PoroAnisotropyStrategyBase> anisotropy_strategy_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
