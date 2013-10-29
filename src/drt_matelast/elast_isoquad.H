/*----------------------------------------------------------------------*/
/*!
\file elast_isoquad.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a quadratic material.

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOQUAD_H_
#define ELAST_ISOQUAD_H_

// included in elast_summand -> elasthyper.H
//#include "../drt_mat/matpar_parameter.H"
// included in elast_summand.H
//#include "../drt_inpar/inpar_material.H"
//#include "../linalg/linalg_fixedsizematrix.H"
#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a cubic material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoQuad C 1
    class IsoQuad : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoQuad(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoQuad() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double c_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetC(const double c) { c_ = c; }
      //@}

    };  // class IsoQuad

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric quadratic material
  ///
  /// This is the isochoric part of a second order hyperelastic, isotropic
  /// material depending on the modified trace of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C (\overline{I}_{\boldsymbol{C}}-3)^2.
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 11/09
  class IsoQuad : public Summand
  {

  public:

    /// empty constructor
    IsoQuad();

    /// constructor with given material parameters
    IsoQuad(MAT::ELASTIC::PAR::IsoQuad* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoquad; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );


    /// @name Access methods
    //@{
    double C() const { return params_->c_; }

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
      isomod = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoQuad* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISOQUAD_H_
