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

      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
          ) = 0;


      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          const double& press,      ///< (i) lubrication pressure at gauss point
          const double& viscosity,  ///< (i) viscosity at gauss point
          double& dviscosity_dp     ///< (o) derivative of viscosity w.r.t. pressure
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

      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
      );

      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          const double& press,      ///< (i) lubrication pressure at gauss point
          const double& viscosity,  ///< (i) viscosity at gauss point
          double& dviscosity_dp     ///< (o) derivative of potential w.r.t. pressure
      );


      /// constant viscosity
      double viscosity_;

    };  // class PoroLawLinear

    /*-------------------------------------------------------------------*/
    //! Lubrication law regarding Barus viscosity
    class LubricationLawBarus : public LubricationLaw
    {
     public:
      /// standard constructor
      explicit LubricationLawBarus(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~LubricationLawBarus() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();


      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
      );

      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          const double& press,      ///< (i) lubrication pressure at gauss point
          const double& viscosity,  ///< (i) viscosity at gauss point
          double& dviscosity_dp     ///< (o) derivative of potential w.r.t. pressure
      );



      /// material parameters for the calculation of the viscosity
      /// Absolute viscosity value in Barus formulation
      double ABSViscosity_;
      /// Press_vis coefficient value in Barus formulation
      double PreVisCoeff_;
    };

    /*-------------------------------------------------------------------*/
    //! Lubrication law regarding Roeland viscosity
    class LubricationLawRoeland : public LubricationLaw
    {
     public:
      /// standard constructor
      explicit LubricationLawRoeland(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~LubricationLawRoeland() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();


      /// compute current viscosity and save it
      virtual void ComputeViscosity(
          const double& press,  ///< (i) lubrication pressure at gauss point
          double& viscosity     ///< (o) viscosity at gauss point
      );

      //! evaluate constitutive relation for viscosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          const double& press,      ///< (i) lubrication pressure at gauss point
          const double& viscosity,  ///< (i) viscosity at gauss point
          double& dviscosity_dp     ///< (o) derivative of potential w.r.t. pressure
      );



      /// material parameters for the calculation of the viscosity
      /// Absolute viscosity value in Roeland formulation
      double ABSViscosity_;
      /// Press_vis coefficient value in Roeland formulation
      double PreVisCoeff_;

      // constants with units for the calculation of the viscosity
      /// Reference viscosity in Roeland formulation
      double RefVisc_;
      /// Reference pressure in Roeland formulation
      double RefPress_;
      /// The analytical term to corrolate Barus and Roeland formulation
      double z_;
    };


  }  // namespace PAR
}  // namespace MAT


#endif /* SRC_DRT_MAT_LUBRICATION_LAW_H_ */
