/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for viscosity for multiphase
 porous flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_FLUIDPORO_VISCOSITY_LAW_HPP
#define FOUR_C_MAT_FLUIDPORO_VISCOSITY_LAW_HPP

#include "4C_config.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    //! interface class for generic viscosity law
    class FluidPoroViscosityLaw : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroViscosityLaw(
          Teuchos::RCP<Core::Mat::PAR::Material> matdata, bool constviscosity)
          : Parameter(matdata), constviscosity_(constviscosity){};

      // get viscosity
      virtual double GetViscosity(const double abspressgrad) const = 0;

      // get derivative of viscosity wrt |grad(p)|
      virtual double get_deriv_of_viscosity_wrt_abs_press_grad(const double abspressgrad) const = 0;

      // check for constant viscosity
      bool has_constant_viscosity() const { return constviscosity_; }

      /// factory method
      static Mat::PAR::FluidPoroViscosityLaw* CreateViscosityLaw(int matID);

     private:
      const bool constviscosity_;
    };

    //! class for constant viscosity law
    class FluidPoroViscosityLawConstant : public FluidPoroViscosityLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroViscosityLawConstant(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; };

      // get viscosity
      double GetViscosity(const double abspressgrad) const override { return viscosity_; };

      // get derivative of viscosity wrt |grad(p)|  --> 0 in case of const. viscosity
      double get_deriv_of_viscosity_wrt_abs_press_grad(const double abspressgrad) const override
      {
        return 0.0;
      };

     private:
      /// @name material parameters
      //@{
      /// viscosity (constant in this case)
      const double viscosity_;
      //@}
    };

    //! class for viscosity law modelling cell adherence from
    // G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and B. A. Schrefler.
    // Three phase flow dynamics in tumor growth. Computational Mechanics, 53:465-484, 2014.
    // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / |
    // grad(pressure) |) + xi)
    class FluidPoroViscosityLawCellAdherence : public FluidPoroViscosityLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroViscosityLawCellAdherence(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; };

      // get viscosity
      double GetViscosity(const double abspressgrad) const override;

      // get derivative of viscosity wrt |grad(p)|
      double get_deriv_of_viscosity_wrt_abs_press_grad(const double abspressgrad) const override;

     private:
      /// @name material parameters
      //@{
      /// viscosity 0
      const double visc0_;
      /// xi
      const double xi_;
      /// psi
      const double psi_;
      //@}
    };


  }  // namespace PAR
}  // namespace Mat



FOUR_C_NAMESPACE_CLOSE

#endif
