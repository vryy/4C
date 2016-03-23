/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke_ActiveStress.H

\brief Active stress material

\maintainer Marc Hirschvogel
*----------------------------------------------------------------------*/


#ifndef ELAST_COUPANISONEOHOOKE_ACTIVESTRESS_H_
#define ELAST_COUPANISONEOHOOKE_ACTIVESTRESS_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    ///
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupAnisoNeoHooke_ActiveStress SIGMA 100.0 TAUC0 0.0 MAX_ACTIVATION 30.0 MIN_ACTIVATION -20.0 SOURCE_ACTIVATION 1 ACTIVATION_THRES 0
    class CoupAnisoNeoHooke_ActiveStress : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAnisoNeoHooke_ActiveStress(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoNeoHooke_ActiveStress() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double sigma_;
      /// initial condition
      double tauc0_;
      /// Maximal value for rescaling the activation curve
      double maxactiv_;
      /// Minimal value for rescaling the activation curve
      double minactiv_;
      /// Threshold for stress activation function
      double activationthreshold_;
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
      void SetC(const double sigma) { sigma_ = sigma; }
      //@}

    };  // class CoupAnisoNeoHooke_ActiveStress

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// This is a simplification of the muscle contraction law proposed in [1], essentially adds an active stress tau
  /// along a given fiber direction f_0
  /// active 2nd Piola-Kirchhoff stress:
  ///\f[
  ///   S_{active} = \tau(t) f_0 \otimes f_0
  ///\f]
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] 2001_Bestel, Clément, Sorine - A Biomechanical Model of Muscle Contraction (2001), Medical Image Computing and
  /// Computer-Assisted Intervention (MICCAI'01), vol. 2208, Springer-Verlag Berlin, 1159-1161
  /// </ul>
  ///
  class CoupAnisoNeoHooke_ActiveStress : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupAnisoNeoHooke_ActiveStress(MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoneohooke_activestress; }

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
    double SIGMA() const { return params_->sigma_; }

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

    // update internal stress variables
    virtual void Update();

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
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress* params_;

    /// fiber direction
    LINALG::Matrix<3,1>  a_;
    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A_;
    /// Active stress at current time step
    double tauc_;
    /// Active stress at last time step
    double tauc_last_;


  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPANISONEOHOOKE_ACTIVESTRESS_H_
