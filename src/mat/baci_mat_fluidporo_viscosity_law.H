/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for viscosity for multiphase
 porous flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_FLUIDPORO_VISCOSITY_LAW_H
#define BACI_MAT_FLUIDPORO_VISCOSITY_LAW_H

#include "baci_config.H"

#include "baci_mat_par_parameter.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic viscosity law
    class FluidPoroViscosityLaw : public Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroViscosityLaw(Teuchos::RCP<MAT::PAR::Material> matdata, bool constviscosity)
          : Parameter(matdata), constviscosity_(constviscosity){};

      // get viscosity
      virtual double GetViscosity(const double abspressgrad) const = 0;

      // get derivative of viscosity wrt |grad(p)|
      virtual double GetDerivOfViscosityWrtAbsPressGrad(const double abspressgrad) const = 0;

      // check for constant viscosity
      bool HasConstantViscosity() const { return constviscosity_; }

      /// factory method
      static MAT::PAR::FluidPoroViscosityLaw* CreateViscosityLaw(int matID);

     private:
      const bool constviscosity_;
    };

    //! class for constant viscosity law
    class FluidPoroViscosityLawConstant : public FluidPoroViscosityLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroViscosityLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      // get viscosity
      double GetViscosity(const double abspressgrad) const override { return viscosity_; };

      // get derivative of viscosity wrt |grad(p)|  --> 0 in case of const. viscosity
      double GetDerivOfViscosityWrtAbsPressGrad(const double abspressgrad) const override
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
      explicit FluidPoroViscosityLawCellAdherence(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      // get viscosity
      double GetViscosity(const double abspressgrad) const override;

      // get derivative of viscosity wrt |grad(p)|
      double GetDerivOfViscosityWrtAbsPressGrad(const double abspressgrad) const override;

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
}  // namespace MAT



BACI_NAMESPACE_CLOSE

#endif  // MAT_FLUIDPORO_VISCOSITY_LAW_H
