/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for (microscopic) density in
 porous media

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_PORO_DENSITY_LAW_H
#define BACI_MAT_PORO_DENSITY_LAW_H

#include "baci_mat_par_parameter.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic density law
    class PoroDensityLaw : public Parameter
    {
     public:
      /// standard constructor
      explicit PoroDensityLaw(Teuchos::RCP<MAT::PAR::Material> matdata) : Parameter(matdata){};

      //! compute derivative of density w.r.t. pressure
      virtual double ComputeCurDensityDerivative(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) = 0;

      /// compute current density
      virtual double ComputeCurDensity(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) = 0;

      /// compute relation of reference density to current density
      virtual double ComputeRefDensityToCurDensity(
          const double& press  ///< (i) pressure at gauss point
          ) = 0;

      /// compute derivative of relation of reference density to current density w.r.t. pressure
      virtual double ComputeRefDensityToCurDensityDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) = 0;
      /// compute second derivative of relation of reference density to current density w.r.t.
      /// pressure
      virtual double ComputeRefDensityToCurDensitySecondDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) = 0;

      /// return inverse bulkmodulus (=compressibility)
      virtual double InvBulkmodulus() const = 0;

      /// factory method
      static MAT::PAR::PoroDensityLaw* CreateDensityLaw(int matID);
    };

    //! class for constant density law
    class PoroDensityLawConstant : public PoroDensityLaw
    {
     public:
      /// standard constructor
      explicit PoroDensityLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata)
          : PoroDensityLaw(matdata){};

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      //! compute derivative of density w.r.t. pressure
      double ComputeCurDensityDerivative(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) override
      {
        return refdensity;
      };

      /// compute current density
      double ComputeCurDensity(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) override
      {
        return 0.0;
      };

      /// compute relation of reference density to current density
      double ComputeRefDensityToCurDensity(const double& press  ///< (i) pressure at gauss point
          ) override
      {
        return 1.0;
      };

      /// compute derivative of relation of reference density to current density w.r.t. pressure
      double ComputeRefDensityToCurDensityDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) override
      {
        return 0.0;
      };
      /// compute second derivative of relation of reference density to current density w.r.t.
      /// pressure
      double ComputeRefDensityToCurDensitySecondDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) override
      {
        return 0.0;
      };

      /// return inverse bulkmodulus (=compressibility)
      double InvBulkmodulus() const override { return 0.0; };
    };

    //! class for exponential density law
    class PoroDensityLawExp : public PoroDensityLaw
    {
     public:
      /// standard constructor
      explicit PoroDensityLawExp(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! compute derivative of density w.r.t. pressure
      double ComputeCurDensityDerivative(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) override;

      /// compute current density
      double ComputeCurDensity(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
          ) override;

      /// compute relation of reference density to current density
      double ComputeRefDensityToCurDensity(const double& press  ///< (i) pressure at gauss point
          ) override;

      /// compute derivative ofrelation of reference density to current density w.r.t. pressure
      double ComputeRefDensityToCurDensityDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) override;
      /// compute second derivative of relation of reference density to current density w.r.t.
      /// pressure
      double ComputeRefDensityToCurDensitySecondDerivative(
          const double& press  ///< (i) pressure at gauss point
          ) override;

      /// return inverse bulkmodulus (=compressibility)
      double InvBulkmodulus() const override { return 1.0 / bulkmodulus_; };

     private:
      /// @name material parameters
      //@{
      /// bulk modulus
      double bulkmodulus_;
      //@}
    };

  }  // namespace PAR
}  // namespace MAT



BACI_NAMESPACE_CLOSE

#endif  // MAT_PORO_DENSITY_LAW_H
