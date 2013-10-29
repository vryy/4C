/*----------------------------------------------------------------------*/
/*!
\file elast_coup1pow.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a general power-type material.

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUP1POW_H_
#define ELAST_COUP1POW_H_

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a general power material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_Coup1Pow C 1 D 1
    class Coup1Pow : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      Coup1Pow(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Coup1Pow() { ; }

      /// @name material parameters
      //@{

      /// material parameters
      double c_;
      int d_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetC(const double c) { c_ = c; }
      void SetD(const int d) { d_ = d; }
      //@}

    };  // class Coup1Pow

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric general power material
  ///
  /// This is a summand of a variable order hyperelastic, isotropic
  /// material depending on the modified trace of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C (I_{\boldsymbol{C}}-3)^D.
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 03/11
  class Coup1Pow : public Summand
  {

  public:

    /// empty constructor
    Coup1Pow();

    /// constructor with given material parameters
    Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coup1pow; }

    //@}

    /// Add principal coeffiencts.
    virtual void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// @name Access methods
    //@{
    double C() const { return params_->c_; }
    int D() const { return params_->d_; }

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
    MAT::ELASTIC::PAR::Coup1Pow* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUP1POW_H_
