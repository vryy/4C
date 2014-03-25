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

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

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

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

      /// @name Write methods
      //@{
      void SetEpsilon(const double eps) { eps_ = eps; }
      void SetGamma(const double gam) { gam_ = gam; }
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

    /// constructor with given material parameters
    VolPenalty(MAT::ELASTIC::PAR::VolPenalty* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_volpenalty; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );

    /// @name Access methods
    //@{
    double Epsilon() const { return params_->eps_; }
    double Gamma() const { return params_->gam_; }

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
   MAT::ELASTIC::PAR::VolPenalty* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // VOLOGDEN_H_
