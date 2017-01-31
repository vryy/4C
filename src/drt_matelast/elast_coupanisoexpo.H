/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpo.H
\brief the input line should read
  MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0

\level 1

<pre>
\maintainer Fabian Braeu
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPO_H_
#define ELAST_COUPANISOEXPO_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"


namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for coupled contribution of a anisotropic exponential fiber material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0
    class CoupAnisoExpo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAnisoExpo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoExpo() { ; }

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

    };  // class CoupAnisoExpo

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Coupled anisotropic exponential fiber function, implemented for one possible fiber family
  /// [1]
  /// This is a hyperelastic, anisotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac {k_1}{2 k_2} \left(e^{k_2 (IV_{\boldsymbol C}-1)^2 }-1 \right).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] G.A. Holzapfel, T.C. Gasser, R.W. Ogden: A new constitutive framework for arterial wall mechanics
  ///          and a comparative study of material models, J. of Elasticity 61 (2000) 1-48.
  /// </ul>
  ///
  /// \author st
  /// \date 03/12

  class CoupAnisoExpo : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupAnisoExpo(MAT::ELASTIC::PAR::CoupAnisoExpo* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoexpo; }

    //@}

    /// Setup of summand
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /// retrieve coefficients of first, second and third derivative
    /// of summand with respect to anisotropic invariants
    ///
    /// The derivatives of the summand \f$\Psi(IV_{\boldsymbol{C},\boldsymbol{a}},V_{\boldsymbol{C},\boldsymbol{a}})\f$
    /// in which the principal anisotropic invariants are the arguments are defined as following:
    ///
    /// First derivatives:
    ///
    ///\f[
    /// dPI_{0,aniso} = \frac{\partial \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    ///\f[
    /// dPI_{1,aniso} = \frac{\partial \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    /// Second derivatives:
    ///\f[
    /// ddPII_{0,aniso} = \frac{\partial^2 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}^2} ;
    ///\f]
    ///\f[
    /// ddPII_{1,aniso} = \frac{\partial^2 \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}^2} ;
    ///\f]
    ///\f[
    /// ddPII_{2,aniso} = \frac{\partial^2 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
    /// \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    /// Third derivatives:
    ///\f[
    /// dddPIII_{0,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
    /// \partial IV_{\boldsymbol{C},\boldsymbol{a}} \partial IV_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    /// dddPIII_{1,aniso} = \frac{\partial^3 \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}
    /// \partial V_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    /// dddPIII_{2,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
    /// \partial IV_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    /// dddPIII_{3,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
    /// \partial V_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
    ///\f]
    template< typename T >
    void GetDerivativesAniso(LINALG::TMatrix<T,2,1>& dPI_aniso,               ///< first derivative with respect to invariants
                             LINALG::TMatrix<T,3,1>& ddPII_aniso,             ///< second derivative with respect to invariants
                             LINALG::TMatrix<T,4,1>& dddPIII_aniso,           ///< third derivative with respect to invariants
                             LINALG::TMatrix<T,3,3> const& rcg,               ///< right Cauchy-Green tensor
                             const int eleGID) const;                         ///< element GID

    /// Add anisotropic principal stresses
    virtual void AddStressAnisoPrincipal(
        const LINALG::Matrix<6,1> rcg,  ///<right Cauchy Green Tensor
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        Teuchos::ParameterList& params, ///< additional parameters for computation of material properties
        const int eleGID                ///< element GID
        );

    /// add strain energy
    virtual void AddStrainEnergy(
      double& psi,                      ///< strain energy functions
      const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
      const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
      const LINALG::Matrix<6,1>& glstrain,///< Green-Lagrange strain in strain like Voigt notation
      const int eleGID                   ///< element GID
      );

    /// Evaluates strain energy for automatic differentiation with FAD
    template< typename T >
    void EvaluateFunc(T& psi,                               ///< strain energy functions
                      LINALG::TMatrix<T,3,3> const& rcg,    ///< Right Cauchy-Green tensor
                      int const eleGID) const;              ///< element GID


    /// @name Access methods
    //@{
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }
    double K1COMP() const { return params_->k1comp_; }
    double K2COMP() const { return params_->k2comp_; }

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
    MAT::ELASTIC::PAR::CoupAnisoExpo* params_;

    /// fiber direction
    LINALG::Matrix<3,1>  a_;
    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_COUPANISOEXPO_H_
