/*----------------------------------------------------------------------*/
/*!
 \brief calculation classes for evaluation of constitutive relation for (mircroscopic) density in
 porous media

   \level 3

   \maintainer  Martin Kronbichler
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_PORO_DENSITY_LAW_H_
#define SRC_DRT_MAT_PORO_DENSITY_LAW_H_

#include "matpar_parameter.H"

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

      /// destructor
      virtual ~PoroDensityLaw() { ; }

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

      /// destructor
      virtual ~PoroDensityLawConstant() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      //! compute derivative of density w.r.t. pressure
      virtual double ComputeCurDensityDerivative(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
      )
      {
        return refdensity;
      };

      /// compute current density
      virtual double ComputeCurDensity(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
      )
      {
        return 0.0;
      };

      /// compute relation of reference density to current density
      virtual double ComputeRefDensityToCurDensity(
          const double& press  ///< (i) pressure at gauss point
      )
      {
        return 1.0;
      };

      /// compute derivative of relation of reference density to current density w.r.t. pressure
      virtual double ComputeRefDensityToCurDensityDerivative(
          const double& press  ///< (i) pressure at gauss point
      )
      {
        return 0.0;
      };
      /// compute second derivative of relation of reference density to current density w.r.t.
      /// pressure
      virtual double ComputeRefDensityToCurDensitySecondDerivative(
          const double& press  ///< (i) pressure at gauss point
      )
      {
        return 0.0;
      };

      /// return inverse bulkmodulus (=compressibility)
      virtual double InvBulkmodulus() const { return 0.0; };
    };

    //! class for exponential density law
    class PoroDensityLawExp : public PoroDensityLaw
    {
     public:
      /// standard constructor
      explicit PoroDensityLawExp(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PoroDensityLawExp() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //! compute derivative of density w.r.t. pressure
      virtual double ComputeCurDensityDerivative(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
      );

      /// compute current density
      virtual double ComputeCurDensity(
          const double& refdensity,  ///< (i) initial/reference density at gauss point
          const double& press        ///< (i) pressure at gauss point
      );

      /// compute relation of reference density to current density
      virtual double ComputeRefDensityToCurDensity(
          const double& press  ///< (i) pressure at gauss point
      );

      /// compute derivative ofrelation of reference density to current density w.r.t. pressure
      virtual double ComputeRefDensityToCurDensityDerivative(
          const double& press  ///< (i) pressure at gauss point
      );
      /// compute second derivative of relation of reference density to current density w.r.t.
      /// pressure
      virtual double ComputeRefDensityToCurDensitySecondDerivative(
          const double& press  ///< (i) pressure at gauss point
      );

      /// return inverse bulkmodulus (=compressibility)
      virtual double InvBulkmodulus() const { return 1.0 / bulkmodulus_; };

     private:
      /// @name material parameters
      //@{
      /// bulk modulus
      double bulkmodulus_;
      //@}
    };

  }  // namespace PAR
}  // namespace MAT



#endif /* SRC_DRT_MAT_PORO_DENSITY_LAW_H_ */
