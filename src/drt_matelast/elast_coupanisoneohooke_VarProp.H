/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke_VarProp.H

\brief This file contains the routines required for the anisotropic material using a neo-Hooke material according to

\level 3

\maintainer Martin Pfaller
            pfaller@lnm.mw.tum.de
            089/289 15264

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISONEOHOOKE_VARPROP_H_
#define ELAST_COUPANISONEOHOOKE_VARPROP_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for anisochoric contribution of a
    /// neo-Hooke material with one fiber direction with space-time varying coefficients
    ///
    ///
    /// <h3>Input line</h3>
    /// MAT 1 CoupAnisoNeoHooke_VarProp C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0
    class CoupAnisoNeoHooke_VarProp : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAnisoNeoHooke_VarProp(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoNeoHooke_VarProp() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double c_;
      /// Where the activation comes from: 0=scatra , >0 Id for FUNCT
      int sourceactiv_;
      /// azimute angle of spherical coordinates
      double gamma_;
      /// polar angle of spherical coordinates
      double theta_;
      /// fiber initalization status
      int init_;
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
      void SetC(const double c) { c_ = c; }
      //@}

    };  // class CoupAnisoNeoHooke_VarProp

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Coupled anisotropic exponential fiber function, implemented for one possible fiber family
  /// [1]
  /// This is a hyperelastic, anisotropic material
  /// of the most simple kind. The material parameter can vary in space and time
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = c(X,t) \left(IV_{\boldsymbol C}-1\right).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, TC Gasser, A viscoelastic model for fiber-reinforced composites at finite strains: ... 2000
  /// </ul>
  ///
  /// \author rausch
  /// \date 05/09
  class CoupAnisoNeoHooke_VarProp : public Summand
  {

  public:

    /// empty constructor
//    CoupAnisoNeoHooke_VarProp();

    /// constructor with given material parameters
    CoupAnisoNeoHooke_VarProp(MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoneohooke_varprop; }

    //@}

    /// Setup of summand
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /// Add anisotropic principal stresses
    virtual void AddStressAnisoPrincipal(
        const LINALG::Matrix<6,1> rcg,  ///<right Cauchy Green Tensor
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        Teuchos::ParameterList& params, ///< additional parameters for computation of material properties
        const int eleGID                ///< element GID
        );
    /// @name Access methods
    //@{
    double C() const { return params_->c_; }

    //@}

    /// Set fiber directions
    virtual void SetFiberVecs(
        const double newgamma,   ///< new angle
        const LINALG::Matrix<3,3> locsys,   ///< local coordinate system
        const LINALG::Matrix<3,3> defgrd   ///< deformation gradient
    );

    /// Get fiber directions
    virtual void GetFiberVecs(
        std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
    );

    /// Setup of patient-specific materials
    virtual void SetupAAA(Teuchos::ParameterList& params, const int eleGID){return;}

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
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp* params_;

    /// fiber direction
    LINALG::Matrix<3,1>  a_;
    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPANISONEOHOOKE_VARPROP_H_
