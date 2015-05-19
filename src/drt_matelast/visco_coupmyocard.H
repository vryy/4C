/*----------------------------------------------------------------------*/
/*!
\file visco_coupmyocard.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a viscous material representing the collagen and elastin matrix surrounding
the myocardial fiber (chappelle12)
<pre>
Maintainer: Martin Pfaller, pfaller@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef VISCO_COUPMYOCARD_H_
#define VISCO_COUPMYOCARD_H_

#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isotropic contribution of myocardial matrix
    ///
    /// <h3>Input line</h3>
    /// MAT 1 VISCO_CoupMyocard N 1
    class CoupMyocard : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupMyocard(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupMyocard() { ; }

      /// @name material parameters
      //@{

      /// material parameters
      double n_;

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
      void SetN(const double n) { n_ = n; }
      //@}

    };  // class CoupMyocard

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric general power material
  ///
  /// This is the isochoric part of a viscohyperelastic, isotropic
  /// material depending on the modified of the rate of the right
  /// Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \eta/2 J_1^2.
  ///\f]
  /// (n = \eta)
  ///
  /// \author pfaller
  /// \date 04/15
  class CoupMyocard : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupmyocard; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsViscoModified(
        const LINALG::Matrix<3,1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
        LINALG::Matrix<8,1>& modmu,  ///<necassary coefficients for piola-kirchhoff-stress
        LINALG::Matrix<33,1>& modxi,  ///<necassary coefficients for viscosity tensor
        LINALG::Matrix<7,1>& modrateinv,
        Teuchos::ParameterList& params,
        const int eleGID
    );

    /// @name Access methods
    //@{
    double N() const { return params_->n_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isovisco, ///< global indicator for isotropic viscous formulation
        bool& viscogenmax  ///< global indicator for viscous contribution according the SLS-Model
      )
    {
      isoprinc = true;
      isovisco = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupMyocard* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPMYOCARD_H_
