/*----------------------------------------------------------------------*/
/*!
\file elast_isomooneyrivlin.H
\brief

This file contains the routines required to calculate the isochoric contribution
of a Mooney-Rivlin-type material

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOMOONEYRIVLIN_H_
#define ELAST_ISOMOONEYRIVLIN_H_

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
    /// material parameters for isochoric contribution of a Mooney-Rivlin type material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoMooneyRivlin C1 100 C2 50
    class IsoMooneyRivlin : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoMooneyRivlin(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoMooneyRivlin() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double c1_;
      double c2_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetC1(const double c1) { c1_ = c1; }
      void SetC2(const double c2) { c2_ = c2; }
      //@}

    };  // class IsoMooneyRivlin

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric Mooney-Rivlin material according to [1]
  ///
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{II}_{\boldsymbol{C}}-3).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk
  /// \date 04/09
  class IsoMooneyRivlin : public Summand
  {

  public:

    /// empty constructor
    IsoMooneyRivlin();

    /// constructor with given material parameters
    IsoMooneyRivlin(MAT::ELASTIC::PAR::IsoMooneyRivlin* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isomooneyrivlin; }

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
    MAT::ELASTIC::PAR::IsoMooneyRivlin* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif
