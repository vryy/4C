/*----------------------------------------------------------------------*/
/*!
\file elast_coupvarga.H
\brief
This file contains the routines required for Varga's material

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef ELAST_COUPVARGA_H_
#define ELAST_COUPVARGA_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


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

    /// empty constructor
    CoupVarga();

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

    /// Main material call to determine  2nd PK stress and material constitutive tensor
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add modified coeffiencts.
    /// Nothing to do for this (coupled) material!
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add anisotropic principal coeffiencts.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add anisotropic modified coeffiencts.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsModifiedAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

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


  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupVarga* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // ELAST_COUPVARGA_H_
