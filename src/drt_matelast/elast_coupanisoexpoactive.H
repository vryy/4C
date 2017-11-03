/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpoactive.H
\brief the input line should read
  MAT 1 ELAST_CoupAnisoExpoActive K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0 S 54000 LAMBDAMAX 1.4 LAMBDA0 0.8 DENS 1050

\level 3

<pre>
\maintainer Fabian Braeu
            braeu@lnm.mw.tum.de
            089/289 15236
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOACTIVE_H_
#define ELAST_COUPANISOEXPOACTIVE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for coupled contribution of an anisotropic active fiber material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupAnisoExpoActive K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0 S 54000 LAMBDAM 1.4 LAMBDA0 0.8 DENS 1050
    class CoupAnisoExpoActive : public MAT::PAR::ParameterAniso
    {
    public:

      /// standard constructor
      CoupAnisoExpoActive(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoExpoActive() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double k1_;
      double k2_;
      /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
      double gamma_;
      /// fiber params for the compressible case
      double k1comp_;
      double k2comp_;
      /// fiber initalization status
      int init_;
      /// adapt angle during remodeling
      int adapt_angle_;
      /// maximum contractile stress
      double s_;
      /// stretch at maximum active force generation
      double lambdamax_;
      /// stretch at zero active force generation
      double lambda0_;
      /// total reference mass density at the beginning of the simulation
      double dens_;

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
      void SetK1(const double k1) { k1_ = k1; }
      void SetK2(const double k2) { k2_ = k2; }
      void SetK1COMP(const double k1comp) { k1comp_ = k1comp; }
      void SetK2COMP(const double k2comp) { k2comp_ = k2comp; }
      //@}

    };  // class CoupAnisoExpoActive

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Coupled anisotropic active fiber function, implemented for one possible fiber family
  ///
  /// This is an active anisotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac {s}{\rho} \left(1.0 + \frac{1}{3} \frac{\left(\lambda_m - 1.0 \right)^3}{\left(\lambda_m - \lambda_0 \right)^2} \right).
  ///\f]
  /// [1]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] Wilson, J.S., S. Baek, and J.D. Humphrey, Parametric study of effects of collagen turnover on the natural history of
  ///          abdominal aortic aneurysms. Proc. R. Soc. A, 2013. 469(2150): p. 20120556.
  /// </ul>
  ///
  /// \author braeu
  /// \date 07/16

  class CoupAnisoExpoActive : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupAnisoExpoActive(MAT::ELASTIC::PAR::CoupAnisoExpoActive* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoexpoactive; }

    //@}

    /// Setup of active summand
    void Setup(DRT::INPUT::LineDefinition* linedef);

    /// retrieve coefficients of first, second and third derivative
    /// of summand with respect to anisotropic invariants
    /// ATTENTION: this is only the passive contribution of the fiber!
    template< typename T >
    void GetDerivativesAniso(LINALG::TMatrix<T,2,1>& dPI_aniso,        ///< first derivative with respect to invariants
                             LINALG::TMatrix<T,3,1>& ddPII_aniso,      ///< second derivative with respect to invariants
                             LINALG::TMatrix<T,4,1>& dddPIII_aniso,    ///< third derivative with respect to invariants
                             LINALG::TMatrix<T,3,3> const& rcg,        ///< right Cauchy-Green tensor
                             const int eleGID) const;                  ///< element GID

    /// Add anisotropic principal stresses
    /// ATTENTION: this is only the passive contribution of the fiber!
    virtual void AddStressAnisoPrincipal(
        const LINALG::Matrix<6,1>& rcg,  ///<right Cauchy Green in "strain-like" Voigt notation
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        Teuchos::ParameterList& params, ///< additional parameters for computation of material properties
        const int eleGID                ///< element GID
        );

    /// Evaluates strain energy for automatic differentiation with FAD
    template< typename T >
    void EvaluateFunc(T& psi,                               ///< strain energy functions
                      LINALG::TMatrix<T,3,3> const& rcg,    ///< Right Cauchy-Green tensor
                      int const eleGID) const;              ///< element GID

    /// evaluate stress and cmat
    /// ATTENTION: this is only the active contribution of the fiber!
    template < typename T >
    void EvaluateActiveStressCmatAniso(LINALG::TMatrix<T,3,3> const& CM,   ///< rigtht Cauchy Green tensor
                                       LINALG::TMatrix<T,6,6> & cmat,      ///< material stiffness matrix
                                       LINALG::TMatrix<T,6,1> & stress,    ///< 2nd PK-stress
                                       const int eleGID) const;            ///< element GID

    // add strain energy
    virtual void AddStrainEnergy(
      double& psi,                      ///< strain energy functions
      const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
      const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
      const LINALG::Matrix<6,1>& glstrain,///< Green-Lagrange strain in strain like Voigt notation
      const int eleGID                   ///< element GID
      );

    /// @name Access methods
    //@{
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }
    double K1COMP() const { return params_->k1comp_; }
    double K2COMP() const { return params_->k2comp_; }

    template < typename T >
    inline void GetDerivativeAnisoActive(T & dPIact) const {
      dPIact = dPIact_;
      return;
    };

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
    MAT::ELASTIC::PAR::CoupAnisoExpoActive* params_;

    /// fiber direction
    LINALG::Matrix<3,1>  a_;

    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A_;

    /// first derivative of active fiber potential w.r.t. the active fiber stretch
    double dPIact_;

    /// active fiber stretch for a given muscle tone
    double lambdaact_;
  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_COUPANISOACTIVEFIBER_H_
