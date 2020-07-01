/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for permeability for multiphase
 porous flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_RELPERMEABILITY_LAW_H_
#define SRC_DRT_MAT_FLUIDPORO_RELPERMEABILITY_LAW_H_

#include "matpar_parameter.H"

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic relative permeability law
    class FluidPoroRelPermeabilityLaw : public Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroRelPermeabilityLaw(
          Teuchos::RCP<MAT::PAR::Material> matdata, bool constrelpermeability)
          : Parameter(matdata), constrelpermeability_(constrelpermeability){};

      /// destructor
      virtual ~FluidPoroRelPermeabilityLaw() { ; }

      // get relative permeability
      virtual double GetRelPermeability(const double saturation) const = 0;

      // get derivative of relative permeability with respect to the saturation of this phase
      virtual double GetDerivOfRelPermeabilityWrtSaturation(const double saturation) const = 0;

      // check for constant permeability
      bool HasConstantRelPermeability() const { return constrelpermeability_; }

      /// factory method
      static MAT::PAR::FluidPoroRelPermeabilityLaw* CreateRelPermeabilityLaw(int matID);

     private:
      const bool constrelpermeability_;
    };

    //! class for constant relative permeability law
    class FluidPoroRelPermeabilityLawConstant : public FluidPoroRelPermeabilityLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroRelPermeabilityLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FluidPoroRelPermeabilityLawConstant() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      // get permeability
      double GetRelPermeability(const double saturation) const { return relpermeability_; };

      // get derivative of relative permeability with respect to the saturation of this phase
      virtual double GetDerivOfRelPermeabilityWrtSaturation(const double saturation) const
      {
        return 0.0;
      };

     private:
      /// @name material parameters
      //@{
      /// permeability (constant in this case)
      const double relpermeability_;
      //@}
    };

    //! class for varying relative permeability
    // relative permeabilty of phase i is calculated via saturation_i^exponent as in
    // G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and B. A. Schrefler.
    // Three phase flow dynamics in tumor growth. Computational Mechanics, 53:465–484, 2014
    class FluidPoroRelPermeabilityLawExponent : public FluidPoroRelPermeabilityLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroRelPermeabilityLawExponent(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FluidPoroRelPermeabilityLawExponent() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      // get permeability
      double GetRelPermeability(const double saturation) const
      {
        if (saturation > minsat_)
          return std::pow(saturation, exp_);
        else
          return std::pow(minsat_, exp_);
      };

      // get derivative of relative permeability with respect to the saturation of this phase
      double GetDerivOfRelPermeabilityWrtSaturation(const double saturation) const
      {
        if (saturation > minsat_)
          return exp_ * std::pow(saturation, (exp_ - 1.0));
        else
          return 0.0;
      };

     private:
      /// @name material parameters
      //@{
      /// exponent
      const double exp_;
      /// minimal saturation used for computation -> this variable can be set by the user to avoid
      /// very small values for the relative permeability
      const double minsat_;
      //@}
    };


  }  // namespace PAR
}  // namespace MAT



#endif /* SRC_DRT_MAT_FLUIDPORO_PERMEABILITY_LAW_H_ */
