/*--------------------------------------------------------------------------*/
/*! \file
\brief calculation classes for evaluation of constitutive relation for lubrication

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_MAT_LUBRICATION_LAW_H_
#define SRC_DRT_MAT_LUBRICATION_LAW_H_


#include "matpar_parameter.H"

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic lubrication law
    class LubricationLaw : public Parameter
    {
     public:
      /// standard constructor
      explicit LubricationLaw(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~LubricationLaw() { ; }

      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          double press,          ///< (i) lubrication pressure at gauss point
          double viscosity,      ///< (i) viscosity at gauss point
          double& dviscosity_dp  ///< (o) derivative of viscosity w.r.t. pressure
          ) = 0;

      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
          ) = 0;
    };

    /*----------------------------------------------------------------------*/
    //! constant lubrication law
    class LubricationLawConstant : public LubricationLaw
    {
     public:
      /// standard constructor
      explicit LubricationLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~LubricationLawConstant() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();


      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          double press,          ///< (i) lubrication pressure at gauss point
          double viscosity,      ///< (i) viscosity at gauss point
          double& dviscosity_dp  ///< (o) derivative of potential w.r.t. pressure
      );

      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
      );

      /// constant viscosity
      double viscosity_;

    };  // class PoroLawLinear


  }  // namespace PAR
}  // namespace MAT


#endif /* SRC_DRT_MAT_LUBRICATION_LAW_H_ */
