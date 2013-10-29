/*----------------------------------------------------------------------*/
/*!
\file elast_isoyeoh.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a Yeoh-type material.

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOYEOH_H_
#define ELAST_ISOYEOH_H_

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
    /// material parameters for isochoric contribution of a Yeoh material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200
    class IsoYeoh : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoYeoh(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoYeoh() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double c1_;
      double c2_;
      double c3_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetC1(const double c1) { c1_ = c1; }
      void SetC2(const double c2) { c2_ = c2; }
      void SetC3(const double c3) { c3_ = c3; }
      //@}

    };  // class IsoYeoh

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric Yeoh material according to [1]
  ///
  /// This is the isochoric part of a third order hyperelastic, isotropic
  /// material depending on the modified trace of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{I}_{\boldsymbol{C}}-3)^2 +
  ///          C3 (\overline{I}_{\boldsymbol{C}}-3)^3.
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch/tk
  /// \date 04/09
  class IsoYeoh : public Summand
  {

  public:

    /// empty constructor
    IsoYeoh();

    /// constructor with given material parameters
    IsoYeoh(MAT::ELASTIC::PAR::IsoYeoh* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoyeoh; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );


    /// @name Access methods
    //@{
    double C1() const { return params_->c1_; }
    double C2() const { return params_->c2_; }
    double C3() const { return params_->c3_; }

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
    MAT::ELASTIC::PAR::IsoYeoh* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISOYEOH_H_
