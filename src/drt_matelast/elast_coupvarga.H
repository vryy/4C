/*----------------------------------------------------------------------*/
/*!
\file elast_coupvarga.H

\brief
This file contains the routines required for Varga's material.
The input line should read
  MAT 1 ELAST_CoupVarga MUE 200.0 BETA 1.0

\level 2

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ELAST_COUPVARGA_H_
#define ELAST_COUPVARGA_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters of Varga's material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupVarga MUE 1.0 BETA 1.0
    class CoupVarga : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupVarga(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupVarga() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double mue_;
      /// 'Anti-modulus'
      double beta_;

      //@}

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };


    };  // class CoupVarga

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Varga's material according to [1], [2]
  ///
  /// This is a compressible, hyperelastic material
  /// of the most simple kind.
  ///
  /// The material strain energy density function is
  ///\f[
  /// \Psi = \underbrace{(2\mu-\beta)}_{\displaystyle\alpha} \Big(\lambda_1 + \lambda_2 + \lambda_3 - 3\Big)
  ///      + \beta \Big(\frac{1}{\lambda_1} + \frac{1}{\lambda_2} + \frac{1}{\lambda_3} - 3\Big)
  ///\f]
  /// which was taken from [1] Eq (6.129) and [2] Eq (1.3).
  ///
  /// The material is stress-free in the \f$\lambda_1=\lambda_2=\lambda_3=1\f$ configuration
  /// if \f$\beta=\alpha=\mu\f$.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// <li> [2] JM Hill and DJ Arrigo, "New families of exact solutions for
  ///          finitely deformed incompressible elastic materials",
  ///          IMA J Appl Math, 54:109-123, 1995.
  /// </ul>
  ///
  /// \author bborn
  /// \date 08/09
  class CoupVarga : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupVarga(MAT::ELASTIC::PAR::CoupVarga* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupvarga; }

    /// add shear modulus equivalent
    virtual void AddShearMod(
      bool& haveshearmod,  ///< non-zero shear modulus was added
      double& shearmod  ///< variable to add upon
      ) const;

    //@}

    /// Answer if coefficients with respect to principal stretches are provided
    virtual bool HaveCoefficientsStretchesPrincipal()
    {
      return true;
    }


    /// Add coefficients with respect to principal stretches (or zeros)
    virtual void AddCoefficientsStretchesPrincipal(
      LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
      LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
      const LINALG::Matrix<3,1>& prstr  ///< principal stretches, [lambda_1, lambda_2, lambda_3]
      );

    /// add the derivatives of a coupled strain energy functions associated with a purely isochoric deformation
    virtual void AddCoupDerivVol(const double j, double& dPj1, double& dPj2, double& dPj3)
    {dserror("not implemented");}

    /// @name Access methods
    //@{
    double Mue() const { return params_->mue_; }
    double Beta() const { return params_->beta_; }
    //@}

    /// @name Write methods
    //@{
    void SetMue(const double& mue) { params_->mue_ = mue; }
    void SetBeta(const double& beta) { params_->beta_ = beta; }
    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
      )
    {
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupVarga* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_COUPVARGA_H_
