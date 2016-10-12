/*----------------------------------------------------------------------*/
/*!
\file elast_anisoactivestress_evolution.H

\brief Active stress material

\maintainer Marc Hirschvogel, originally by C. Bertoglio

\level 2
*----------------------------------------------------------------------*/


#ifndef ELAST_ANISOACTIVESTRESS_EVOLUTION_H_
#define ELAST_ANISOACTIVESTRESS_EVOLUTION_H_

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
    /// MAT 1 ELAST_AnisoActiveStress_Evolution SIGMA 100.0 TAUC0 0.0 MAX_ACTIVATION 30.0 MIN_ACTIVATION -20.0 SOURCE_ACTIVATION 1 ACTIVATION_THRES 0 [STRAIN_DEPENDENCY No] [LAMBDA_LOWER 0.707] [LAMBDA_UPPER 1.414]
    class AnisoActiveStress_Evolution : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      AnisoActiveStress_Evolution(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~AnisoActiveStress_Evolution() { ; }

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
      /// is there a strain dependency for the active tension?
      int strain_dep_;
      /// lower stretch threshold
      double lambda_lower_;
      /// upper stretch threshold
      double lambda_upper_;
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

    };  // class AnisoActiveStress_Evolution

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// This is a simplification of the muscle contraction law proposed in [1],[2],
  /// resulting in the following first order ODE for the active stress tau, compare [3]:
  ///\f[
  ///   \frac{d}{dt} \tau = n_0 \sigma_0 |u|_+ - \tau |u|, \quad \tau(0) = tau_0
  ///\f]
  /// where \f$\sigma_0\f$ is the contractility (asymptotic value of \tau) and u is a control variable either provided by
  /// a electrophysiology simulation or by a user-specified function
  /// therein, n0 is a strain-dependent factor that may take into account the Frank-Starling effect!
  /// n0 \in [0; 1] scales the contractility depending on a lower and upper fiber stretch lambda
  /// using a flipped parabola, n0 = -(lambda - lambda_lower)*(lambda - lambda_upper) * 4/(lambda_lower-lambda_upper)^2,
  /// which is an approximation to the function from Sainte-Marie et al. 2006, Fig. 2(ii)
  /// other laws might be thought of here, since hardly any literature provides a meaningful dependency...
  ///
  /// Due to the active stress approach, see [4], the active stress will be added along a given
  /// fiber direction f_0 to the 2nd Piola-Kirchhoff stress:
  ///\f[
  ///   S_{active} = \tau(t) f_0 \otimes f_0
  ///\f]
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] 2012 Chapelle, Le Tallec, Moireau, Sorine - An energy-preserving muscle tissue model: formulation and compatible discretizations,
  // Journal for Multiscale Computational Engineering 10(2):189-211 (2012)
  /// <li> [2] 2001 Bestel, Clément, Sorine - A Biomechanical Model of Muscle Contraction (2001), Medical Image Computing and
  /// Computer-Assisted Intervention (MICCAI'01), vol. 2208, Springer-Verlag Berlin, 1159-1161
  /// <li> [3] 2002 Sermesant, Coudièr, Delingette, Ayache - Progress towards an electromechanical model of the heart for cardiac image analysis. (2002)
  /// IEEE International Symposium on Biomedical Imaging, 10-13
  /// <li> [4] 1998 Hunter, McCulloch, ter Keurs - Modelling the mechanical properties of cardiac muscle (1998),
  /// Progress in Biophysics and Molecular Biology
  /// </ul>
  ///
  class AnisoActiveStress_Evolution : public Summand
  {

  public:

    /// constructor with given material parameters
    AnisoActiveStress_Evolution(MAT::ELASTIC::PAR::AnisoActiveStress_Evolution* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_anisoactivestress_evolution; }

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
    MAT::ELASTIC::PAR::AnisoActiveStress_Evolution* params_;

    /// fiber direction
    LINALG::Matrix<3,1>  a_;
    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A_;
    /// Active stress at current time step
    double tauc_np_;
    /// Active stress at last time step
    double tauc_n_;


  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ANISOACTIVESTRESS_EVOLUTION_H_
