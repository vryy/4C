/*----------------------------------------------------------------------*/
/*!
\file elast_volpenalty.H
\brief
This file contains the routines required for the volumetic function
as afore used in the lungpenaly material for more details see
[Roernbauer2008]

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLPENALTY_H_
#define ELAST_VOLPENALTY_H_
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
    /// material parameters for volumetric contribution
    /// \f$\Psi=\epsilon \left( J^{\gamma} + \frac 1 {J^{\gamma}} -2 \right)\f$
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_VolPenalty EPSILON 1. GAMMA 1.
    class VolPenalty : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VolPenalty(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VolPenalty() { ; }

      /// @name material parameters
      //@{

      /// Dilatation modulus
      double eps_;
      /// empiric constant
      double gam_;

      //@}

    };  // class VolPenalty

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Volumetric Penalty
  ///
  /// Strain energy function is given by
  ///
  ///
  ///\f[
  ///   \Psi=\epsilon \left( J^{\gamma} + \frac 1 {J^{\gamma}} -2 \right)
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 08/09
  class VolPenalty : public Summand
  {

  public:

    /// empty constructor
    VolPenalty();

    /// constructor with given material parameters
    VolPenalty(MAT::ELASTIC::PAR::VolPenalty* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_volpenalty; }

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

        /// @name Access methods
    //@{
    double Epsilon() const { return params_->eps_; }
    double Gamma() const { return params_->gam_; }

    //@}

    /// @name Write methods
    //@{
    void SetEpsilon(const double eps) { params_->eps_ = eps; }
    void SetGamma(const double gam) { params_->gam_ = gam; }
    //@}


  private:

    /// my material parameters
   MAT::ELASTIC::PAR::VolPenalty* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // VOLOGDEN_H_
