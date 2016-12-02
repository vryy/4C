/*----------------------------------------------------------------------*/
/*!
\file elast_coupsimopister.H

\brief
This file contains the routines required to calculate the Simo and Pister
material model. (U(J) is not implemented).
Always use in combination with at least one other material model of the
Elasthyper Toolbox.
The input line should read
  MAT 1 ELAST_CoupSimoPister MUE 1000

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPSIMOPISTER_H_
#define ELAST_COUPSIMOPISTER_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a CoupSimoPisteran material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupSimoPister MUE 1000
    class CoupSimoPister : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupSimoPister(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupSimoPister() { ; }

      /// @name material parameters
      //@{

      /// material parameters


      /// Material Parameter my
      double mue_;


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
      void SetMue(const double mue) { mue_ = mue; }
      //@}

    };  // class CoupSimoPister

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// CoupSimoPisteran material
  ///
  /// This is the summand of a hyperelastic, isotropic CoupSimoPisteran material
  /// depending on the first and the third invariant of the right Cauchy-Green tensor.
  /// The formulation is based on [1] (4.3) or [2] (1)
  ///
  /// The implemented material is a coupled form of the Simo & Pister material
  /// model. The part U(J) is not implemented. So this material model can/should
  /// combined with any volumetric material model.
  /// The Parameter read in is my.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = 0.5*\mu(I_{\boldsymbol{C}}-3) - \mu log(J)
  ///\f]
  ///
  /// [1] Simo and Pister - 1984
  /// [2] Hartmann - "The class of Simo & Pister-type Hyperelasticity Relastions"
  ///
  /// \author birzle
  /// \date 04/16
  class CoupSimoPister : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupSimoPister(MAT::ELASTIC::PAR::CoupSimoPister* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupsimopister; }

    //@}

    // add strain energy
    virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy function
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const int eleGID                   ///< element GID
        );

    virtual void AddDerivativesPrincipal(
        LINALG::Matrix<3,1>& dPI, ///< first derivative with respect to invariants
        LINALG::Matrix<6,1>& ddPII, ///< second derivative with respect to invariants
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const int eleGID                  ///< element GID
        );

    /// @name Access methods
    //@{
    double Mue() const { return params_->mue_; }

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
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupSimoPister* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPSIMOPISTER_H_
