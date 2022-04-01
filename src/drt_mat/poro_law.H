/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for porosity


\level 2
 *----------------------------------------------------------------------*/


#ifndef PORO_LAW_H_
#define PORO_LAW_H_

#include "matpar_parameter.H"

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic porosity law
    class PoroLaw : public Parameter
    {
     public:
      //! standard constructor
      explicit PoroLaw(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! evaluate constitutive relation for porosity and compute derivatives
      virtual void ConstitutiveDerivatives(
          const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,                   //!< (i) fluid pressure at gauss point
          const double& J,                       //!< (i) Jacobian determinant at gauss point
          const double& porosity,                //!< (i) porosity at gauss point
          const double& refporosity,             //!< (i) porosity at gauss point
          double* dW_dp,                         //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,                       //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,                         //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,  //!< (o) derivative of potential w.r.t. reference porosity
          double* W            //!< (o) inner potential
          ) = 0;

      //! compute current porosity and save it
      virtual void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) = 0;

      //! return inverse bulkmodulus (=compressibility)
      virtual double InvBulkModulus() const = 0;
    };

    /*----------------------------------------------------------------------*/
    //! linear porosity law
    class PoroLawLinear : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawLinear(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;


      //! evaluate constitutive relation for porosity and compute derivatives using reference
      //! porosity
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override { return 1.0 / bulk_modulus_; }

     private:
      //! @name material parameters
      //!@{
      //! bulk modulus of skeleton phase
      double bulk_modulus_;
      //!@}
    };

    /*----------------------------------------------------------------------*/
    //! Neo-Hookean like porosity law
    // see   A.-T. Vuong, L. Yoshihara, W.A. Wall :A general approach for modeling interacting flow
    // through porous media under finite deformations, Comput. Methods Appl. Mech. Engrg. 283 (2015)
    // 1240–1259, equation (39)

    class PoroLawNeoHooke : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! evaluate constitutive relation for porosity and compute derivatives
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override { return 1.0 / bulk_modulus_; }

     private:
      //! @name material parameters
      //!@{
      //! bulk modulus of skeleton phase
      double bulk_modulus_;
      //! penalty parameter for porosity
      double penalty_parameter_;
      //!@}
    };

    /*----------------------------------------------------------------------*/
    //! constant porosity law
    class PoroLawConstant : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! evaluate constitutive relation for porosity and compute derivatives
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override { return 0.0; }
    };

    /*----------------------------------------------------------------------*/
    //! incompressible skeleton porosity law
    class PoroLawIncompSkeleton : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawIncompSkeleton(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! evaluate constitutive relation for porosity and compute derivatives
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override { return 0.0; }
    };

    /*----------------------------------------------------------------------*/
    //! linear Biot model for porosity law
    class PoroLawLinBiot : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawLinBiot(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! evaluate constitutive relation for porosity and compute derivatives
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override { return inv_biot_modulus_; }

     private:
      //! @name material parameters
      //!@{
      //! inverse biot modulus
      double inv_biot_modulus_;
      //! Biot coefficient
      double biot_coeff_;
      //!@}
    };

    class PoroDensityLaw;
    /*----------------------------------------------------------------------*/
    //! porosity law depending on density
    class PoroLawDensityDependent : public PoroLaw
    {
     public:
      //! standard constructor
      explicit PoroLawDensityDependent(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! evaluate constitutive relation for porosity and compute derivatives
      void ConstitutiveDerivatives(const Teuchos::ParameterList& params,  //!< (i) parameter list
          const double& press,        //!< (i) fluid pressure at gauss point
          const double& J,            //!< (i) Jacobian determinant at gauss point
          const double& porosity,     //!< (i) porosity at gauss point
          const double& refporosity,  //!< (i) porosity at gauss point
          double* dW_dp,              //!< (o) derivative of potential w.r.t. pressure
          double* dW_dphi,            //!< (o) derivative of potential w.r.t. porosity
          double* dW_dJ,              //!< (o) derivative of potential w.r.t. jacobian
          double* dW_dphiref,         //!< (o) derivative of potential w.r.t. reference porosity
          double* W                   //!< (o) inner potential
          ) override;

      //! compute current porosity and save it
      void ComputePorosity(
          const double& refporosity,  //!< (i) initial/reference porosity at gauss point
          const double& press,        //!< (i) pressure at gauss point
          const double& J,            //!< (i) determinant of jacobian at gauss point
          const int& gp,              //!< (i) number of current gauss point
          double& porosity,           //!< (o) porosity at gauss point
          double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss
                              //!< point
          double* dphi_dJJ,   //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
          double* dphi_dpp,   //!< (o) second derivative of porosity w.r.t. pressure at gauss point
          double* dphi_dphiref  //!< (o) derivative of porosity w.r.t. reference porosity (only
                                //!< nonzero with reaction)
          ) override;

      //! return inverse bulkmodulus (=compressibility)
      double InvBulkModulus() const override;

     private:
      //! @name material parameters
      //!@{
      //! density law
      MAT::PAR::PoroDensityLaw* density_law_;
      //!@}
    };
  }  // namespace PAR
}  // namespace MAT

#endif /* PORO_LAW_H_ */
