/*----------------------------------------------------------------------*/
/*!
\file elast_coupblatzko.H
\brief
This file contains the routines required for Blatz and Ko material model
according to Holzapfel, "Nonlinear solid mechanics", 2001.

<pre>
Maintainer: Sophie Rausch
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPBLATZKO_H_
#define ELAST_COUPBLATZKO_H_

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
    /// material parameters for Blatz and Ko material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupBlatzKo MUE 1.044E7 NUE 0.3 F 0.5
    class CoupBlatzKo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupBlatzKo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupBlatzKo() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double mue_;
      /// Possion's ratio
      double nue_;
      /// interpolation parameter
      double f_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetMue(const double mue) { mue_ = mue; }
      void SetNue(const double nue) { nue_ = nue; }
      void SetF(const double f) { f_ = f; }
      //@}

    };  // class CoupBlatzKo

  }  // namespace PAR
  /*----------------------------------------------------------------------*/
  /// Coupled Blatz and Ko material according to [1].
  ///
  /// Strain energy function is given by
  ///
  ///
  ///\f[
  ///   \Psi= f \frac {\mu} 2 \left[ (I_{\boldsymbol C}-3)+\frac 1
  ///         {\beta} ( III_{\boldsymbol C}^{-\beta} -1) \right]
  ///         +(1-f) \frac {\mu} 2 \left[\left( \frac {II_{\boldsymbol
  ///         C}}{III_{\boldsymbol C}}-3 \right) + \frac 1 {\beta}
  ///         (III_{\boldsymbol C}^{\beta}-1)\right]
  ///\f]
  ///
  /// \f[
  ///    with \ \beta=\nu-0.5
  /// \f]
  ///
  ///
  /// [1] Holzapfel, G.A. "Nonlinear Solid Mechanics", 2000, p.247
  /// \author rausch
  /// \date 04/09
  class CoupBlatzKo : public Summand
  {

  public:

    /// empty constructor
    CoupBlatzKo();

    /// constructor with given material parameters
    CoupBlatzKo(MAT::ELASTIC::PAR::CoupBlatzKo* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupblatzko; }
    //@}

    /// Add principal coeffiencts.
    virtual void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// @name Access methods
    //@{
    double Mue() const { return params_->mue_; }
    double Nue() const { return params_->nue_; }
    double F() const { return params_->f_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco ///< global indicator for isotropic splitted viscos formulation
      )
    {
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupBlatzKo* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPBLATZKO_H_
