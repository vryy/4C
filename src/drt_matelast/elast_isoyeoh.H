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
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

namespace MAT
{
namespace ELAST
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
      const double c1_;
      const double c2_;
      const double c3_;

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
    IsoYeoh(MAT::ELAST::PAR::IsoYeoh* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoyeoh; }

    //@}

    /// Add principal coeffiencts.
    /// Nothing to do for this (decoupled) material!
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
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );

  private:

    /// my material parameters
    MAT::ELAST::PAR::IsoYeoh* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // ISOYEOH_H_
