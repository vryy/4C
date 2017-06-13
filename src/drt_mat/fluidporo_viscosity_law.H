/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_viscosity_law.H

 \brief calculation classes for evaluation of constitutive relation for viscosity for multiphase porous flow

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_VISCOSITY_LAW_H_
#define SRC_DRT_MAT_FLUIDPORO_VISCOSITY_LAW_H_

#include "matpar_parameter.H"

namespace MAT
{
namespace PAR
{
  //! interface class for generic viscosity law
  class FluidPoroViscosityLaw : public Parameter
  {

  public:

    /// standard constructor
    explicit FluidPoroViscosityLaw(
      Teuchos::RCP<MAT::PAR::Material> matdata,
      bool constviscosity
      ):Parameter(matdata),
      constviscosity_(constviscosity)
    {};

    /// destructor
    virtual ~FluidPoroViscosityLaw() { ; }

    // get viscosity
    virtual double GetViscosity(const double abspressgrad) const=0;

    // check for constant viscosity
    bool HasConstantViscosity() const
    {return constviscosity_;}

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
     explicit FluidPoroViscosityLawConstant(
       Teuchos::RCP<MAT::PAR::Material> matdata
       );

     /// destructor
     virtual ~FluidPoroViscosityLawConstant() { ; }

     /// create material instance of matching type with my parameters
     virtual Teuchos::RCP<MAT::Material> CreateMaterial()
            {return Teuchos::null;};

     // get viscosity
     double GetViscosity(const double abspressgrad) const
     {return viscosity_;};

   private:
    /// @name material parameters
    //@{
    /// viscosity (constant in this case)
    const double viscosity_;
    //@}

   };

   //! class for viscosity law modelling cell adherence from
   // G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and B. A. Schrefler.
   // Three phase flow dynamics in tumor growth. Computational Mechanics, 53:465–484, 2014.
   // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / | grad(pressure) |) + xi)
    class FluidPoroViscosityLawCellAdherence : public FluidPoroViscosityLaw
    {

    public:

      /// standard constructor
      explicit FluidPoroViscosityLawCellAdherence(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~FluidPoroViscosityLawCellAdherence() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial()
             {return Teuchos::null;};

      // get viscosity
      double GetViscosity(const double abspressgrad) const;

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


}
}



#endif /* SRC_DRT_MAT_FLUIDPORO_VISCOSITY_LAW_H_ */
