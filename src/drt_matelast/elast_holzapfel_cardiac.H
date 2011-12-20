/*----------------------------------------------------------------------*/
/*!
\file elast_holzapfel_cardiac.H
\brief
This file contains the routines required for the passive behaviour of cardiac muscle
according to Holzapfel and Ogden, "Constitutive modelling of passive myocardium", 2009.

<pre>
Maintainer: Andreas Nagler
            nagler@lnm.mw.tum.de
            089/289 15249
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_HOLZAPFEL_CARDIAC_H_
#define ELAST_HOLZAPFEL_CARDIAC_H_
#ifdef CCADISCRET

// Linalg matrices used instead; included in elast_summand
//#include <Epetra_SerialDenseVector.h>
//#include <Epetra_SerialDenseMatrix.h>
// included in elast_summand
//#include "../drt_inpar/inpar_material.H"
#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a passive cardiac material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_Holzapfel_Cardiac A 57 B 8.023 A4 18472 B4 16.026 A6 2.481 B6 11.120 A8 216 B8 11.436
    class Holzapfel_Cardiac : public MAT::PAR::Parameter
    {
      public:

      /// constructor with given material parameters
      Holzapfel_Cardiac(
              Teuchos::RCP<MAT::PAR::Material> matdata
              );

      /// destructor
      virtual ~Holzapfel_Cardiac() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      const double A_;
      const double B_;
      const double A4_;
      const double B4_;
      const double A6_;
      const double B6_;
      const double A8_;
      const double B8_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      //void SetA(const double A) { A_ = A; }
      //void SetB(const double B) { B_ = B; }
      //void SetK3(const double k3) { k3_ = k3; }
      //void SetK4(const double k4) { k4_ = k4; }
      //@}

    };  // class HolzapfelCardiac

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Anisotropic cardiac material, implemented with two possible fiber families acc [1]
  ///
  /// This is a hyperelastic, anisotropic material for the passive response of cardiac material
  ///
  /// Strain energy function is given by:
  ///\f[
  ///   \Psi = \frac {a}{2 b} \left( e^{b (I_{\boldsymbol C} - 3)} \right) + \frac {a_4}{2 b_4} \left( e^{b_4 (IV_{\boldsymbol C} - 1)^2} - 1 \right) +
  ///   \frac {a_6}{2 b_6} \left( e^{b_6 (VI_{\boldsymbol C} - 1)^2} - 1 \right) + \frac {a_8}{2 b_8} \left( e^{b_8 VIII_{\boldsymbol C}^2} - 1 \right)
  ///\f]
  ///
  /// The isotropic and anisotropic part of the strain energy function \Psi are splitted and treated separately as proposed in [2].
  /// The isotropic part:
  ///\f[
  ///     \Psi_{iso} = \frac {a}{2 b} \left(e^{b (I_{\boldsymbol C} - 3)} \right)
  /// \f]
  /// will be solved in the decoupled variant,i.e.
  ///
  ///\f[
  ///     \overline{\Psi}_{iso} = \frac {a}{2 b} \left(e^{b (J^{\frac{2}{3}} \overline{I}_{\boldsymbol C} - 3)} \right)
  ///\f]
  ///
  /// The anisotropic part:
  ///\f[
  ///   \Psi_{aniso} = \frac {a_4}{2 b_4} \left( e^{b_4 (IV_{\boldsymbol C} - 1)^2} - 1 \right) + \frac {a_6}{2 b_6} \left( e^{b_6 (VI_{\boldsymbol C} - 1)^2} - 1 \right) + \frac {a_8}{2 b_8} \left( e^{b_8 VIII_{\boldsymbol C}^2} - 1 \right)
  ///\f]
  /// will be calculated normal
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, RW Ogden, Constitutive modelling of passive myocardium: a structurally based framework for material characterization
  /// </ul>
  /// <ul>
  /// <li> [2] C Sansour, On the physical assumptions underlying the volumatric-isochoric split and the case of anisotropy
  /// </ul>
  ///
  /// \author nagler
  /// \date 11/11
  class Holzapfel_Cardiac : public Summand
  {

  public:

    /// empty constructor
    Holzapfel_Cardiac();

    /// constructor with given material parameters
    Holzapfel_Cardiac(MAT::ELASTIC::PAR::Holzapfel_Cardiac* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_holzapfel_cardiac; }

    //@}

    /// Add principal coeffiencts.
    /// Nothing to do for this  material!
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;

      return;
    };

    /// Add modified coeffiencts.
    /// Isotropic part will be implemented via the modified implementation
    ///
    /// Coefficients of Piola-Kirchhoff-stress
    ///\f[
    ///   \overline{\Psi}_{iso,1} = \frac{\partial \overline{\Psi}_{iso}}{\partial \overline{I}_{\boldsymbol C}} = a J^{\frac{2}{3}} e^{b (J^{\frac{2}{3}} \overline{I}_{\boldsymbol C} - 3)}
    ///\f]
    ///\f[
    ///   \overline{\Psi}_{iso,2} = \frac{\partial \overline{\Psi}_{iso}}{\partial \overline{II}_{\boldsymbol C}} = 0
    ///\f]
    ///\f[
    ///   \overline{\gamma}_1 = 2  \overline{\Psi}_{iso,1} + 2  \overline{I}_{\boldsymbol C}  \overline{\Psi}_{iso,2} = 2  \overline{\Psi}_{iso,1}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_2 = 0
    ///\f]
    ///
    /// Coefficients of elasticity tensor
    ///\f[
    ///   \overline{\Psi}_{iso,11} = \frac{\partial^2 \overline{\Psi}_{iso}}{\partial \overline{I}_{\boldsymbol C}^2} = 0.5 a b J^{\frac{4}{3}} e^{b (J^{\frac{2}{3}} \overline{I}_{\boldsymbol C} - 3)}
    ///\f]
    ///\f[
    ///   \overline{\Psi}_{iso,12} = \overline{\Psi}_{iso,21} = \overline{\Psi}_{iso,22} = 0
    ///\f]
    ///\f[
    ///   \overline{\gamma}_1 = 4  \left(\overline{\Psi}_{iso,11} + 2 \overline{I}_{\boldsymbol C} \overline{\Psi}_{iso,12} + \overline{I}_{\boldsymbol C}^2 \overline{\Psi}_{iso,22} + \overline{\Psi}_{iso,2} \right) =  2 a b J^{\frac{4}{3}} e^{b (J^{\frac{2}{3}} \overline{I}_{\boldsymbol C} - 3)}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_2 = \overline{\gamma}_3 = \overline{\gamma}_4 = 0
    ///\f]
    ///
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );

    /// Add anisotropic principal coeffiencts.
    /// Anisotropic part will be implemented normal
    ///
    /// Coefficients of Piola-Kirchhoff-stress
    ///\f[
    ///   \Psi_{aniso,4} = a_4 (IV_{\boldsymbol C} - 1) e^{b_4 (IV_{\boldsymbol C} - 1)^2}
    ///\f]
    ///\f[
    ///   \Psi_{aniso,6} = a_6 (VI_{\boldsymbol C} - 1) e^{b_6 (VI_{\boldsymbol C} - 1)^2}
    ///\f]
    ///\f[
    ///   \Psi_{aniso,8} = a_8 VIII_{\boldsymbol C} e^{b_8 VIII_{\boldsymbol C} ^2}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_1^A = 2 \Psi_{aniso,4}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_2^A = 2 \Psi_{aniso,6}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_3^A = \Psi_{aniso,8}
    ///\f]
    ///
    /// Coefficients of elasticity tensor
    ///\f[
    ///   \Psi_{aniso,44} = a_4 e^{b_4 (IV_{\boldsymbol C} - 1)^2} \left( 1 + 2 b_4 (IV_{\boldsymbol C} - 1)^2 \right)
    ///\f]
    ///\f[
    ///   \Psi_{aniso,66} = a_6 e^{b_6 (VI_{\boldsymbol C} - 1)^2} \left( 1 + 2 b_6 (VI_{\boldsymbol C} - 1)^2 \right)
    ///\f]
    ///\f[
    ///   \Psi_{aniso,88} = a_8 e^{b_8 VIII_{\boldsymbol C}^2} \left( 1 + 2 b_8 VIII_{\boldsymbol C}^2 \right)
    ///\f]
    ///\f[
    ///   \overline{\gamma}_1^A = 4 \Psi_{aniso,44}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_2^A = 4 \Psi_{aniso,66}
    ///\f]
    ///\f[
    ///   \overline{\gamma}_{15}^A = \Psi_{aniso,88}
    ///\f]
    ///
    void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    );

    /// @name Access methods
    //@{
    double A() const { return params_->A_; }
    double B() const { return params_->B_; }
    //double K3() const { return params_->k3_; }
    //double K4() const { return params_->k4_; }
    //@}

    private:

    /// my material parameters
    MAT::ELASTIC::PAR::Holzapfel_Cardiac* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // HOLZAPFEL_CARDIAC_H_
