/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpotwocoup.H

\brief This file contains the routines required for the passive behaviour of cardiac muscle according to Holzapfel and Ogden, "Constitutive modelling of passive myocardium", 2009.

\maintainer Marc Hirschvogel, originally by A. Nagler

\level 2

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOTWOCOUP_H_
#define ELAST_COUPANISOEXPOTWOCOUP_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for coupled passive cardiac material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481 B6 11.120 A8 0.216 B8 11.436 GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]
    class CoupAnisoExpoTwoCoup : public MAT::PAR::ParameterAniso
    {
      public:

      /// constructor with given material parameters
      CoupAnisoExpoTwoCoup(
              Teuchos::RCP<MAT::PAR::Material> matdata
              );

      /// destructor
      virtual ~CoupAnisoExpoTwoCoup() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double A4_;
      double B4_;
      double A6_;
      double B6_;
      double A8_;
      double B8_;
      /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
      double gamma_;
      /// fiber initalization status
      int init_;
      /// fibers support compression - or not
      int fib_comp_;
      /// adapt angle during remodeling
      int adapt_angle_;

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
      void SetA4(const double A4) { A4_ = A4; }
      void SetB4(const double B4) { B4_ = B4; }
      void SetA6(const double A6) { A6_ = A6; }
      void SetB6(const double B6) { B6_ = B6; }
      void SetA8(const double A8) { A8_ = A8; }
      void SetB8(const double B8) { B8_ = B8; }
      //@}

    };  // class CoupAnisoExpoTwoCoup

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Anisotropic cardiac material, implemented with two possible fiber families acc [1]
  ///
  /// This is a hyperelastic, anisotropic material for the passive response of cardiac material
  ///
  /// Strain energy function is given by:
  ///\f[
  ///   \Psi = \frac {a_4}{2 b_4} \left( e^{b_4 (IV_{\boldsymbol C} - 1)^2} - 1 \right) +
  ///   \frac {a_6}{2 b_6} \left( e^{b_6 (VI_{\boldsymbol C} - 1)^2} - 1 \right) + \frac {a_8}{2 b_8} \left( e^{b_8 \left( VIII_{\boldsymbol C} - a_0 \cdot b_0 \right)^2} - 1 \right)
  ///\f]
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, RW Ogden, Constitutive modelling of passive myocardium: a structurally based framework for material characterization
  /// </ul>#include "../drt_mat/matpar_parameter.H"
  /// <ul>
  /// <li> [2] C Sansour, On the physical assumptions underlying the volumetric-isochoric split and the case of anisotropy
  /// </ul>
  ///
  /// \author nagler
  /// \date 11/11
  class CoupAnisoExpoTwoCoup : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupAnisoExpoTwoCoup(MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoexpotwocoup; }

    //@}

    /// Setup of summand
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /// Add anisotropic principal stresses
    virtual void AddStressAnisoPrincipal(
        const LINALG::Matrix<6,1>& rcg,  ///<right Cauchy Green Tensor
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        Teuchos::ParameterList& params, ///< additional parameters for computation of material properties
        const int eleGID                ///< element GID
        );
    /// @name Access methods
    //@{
    double A4() const { return params_->A4_; }
    double B4() const { return params_->B4_; }
    double A6() const { return params_->A6_; }
    double B6() const { return params_->B6_; }
    double A8() const { return params_->A8_; }
    double B8() const { return params_->B8_; }
    //@}

    /// Set fiber directions
    virtual void SetFiberVecs(
        const double newgamma,   ///< new angle
        const LINALG::Matrix<3,3>& locsys,   ///< local coordinate system
        const LINALG::Matrix<3,3>& defgrd   ///< deformation gradient
    );

    /// Get fiber directions
    virtual void GetFiberVecs(
        std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
    );

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
      )
    {
      anisoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params_;

    /// fiber directions
    LINALG::Matrix<3,1>  a1_;
    LINALG::Matrix<3,1>  a2_;
    /// dot product fiber direction
    double a1a2_;
    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A1_;
    LINALG::Matrix<6,1>  A2_;
    /// mixed structural tensor (asymmetric) \f$a1 \otimes a2\f$
    LINALG::Matrix<3,3>  A1A2_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPANISOEXPOTWOCOUP_H_
