/*----------------------------------------------------------------------*/
/*! \file
\brief Weickenmeier active skeletal muscle material model

\level 3
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  definitions                                                         |
 *----------------------------------------------------------------------*/
#ifndef MUSCLE_WEICKENMEIER_H
#define MUSCLE_WEICKENMEIER_H

/*----------------------------------------------------------------------*
 |  Headers                                                             |
 *----------------------------------------------------------------------*/
#include "../headers/FADmatrix_utils.H"
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/anisotropy_extension_default.H"
#include "../drt_mat/anisotropy_extension_provider.H"

/*---------------------------------------------------------------------------*
 | Class definitions                                                         |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  // forward declaration
  class Muscle_Weickenmeier;

  namespace PAR
  {
    class Muscle_Weickenmeier : public Parameter
    {
      friend class MAT::Muscle_Weickenmeier;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      Muscle_Weickenmeier(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// passive material parameters
      const double alpha_;   /// material parameter, >0
      const double beta_;    /// material parameter, >0
      const double gamma_;   /// material parameter, >0
      const double kappa_;   /// material parameter for coupled volumetric contribution
      const double omega0_;  /// weighting factor for isotropic tissue constituents, governs ratio
                             /// between muscle matrix material (omega0) and muscle fibers (omegap)
                             /// with omega0 + omegap = 1

      /// active microstructural parameters
      // stimulation frequency dependent activation contribution
      const double
          Na_;  /// number of active motor units (MU) per undeformed muscle cross-sectional area
      const int muTypesNum_;  /// number of motor unit (MU) types
                              /// vectors indicating corresponding parameters of motor unit types
                              /// e.g. slow type I fibres (index 0), fast resistant  type IIa fibres
                              /// (index 1), fast fatigue type IIb fibres (index 2)
      const std::vector<double> I_;    /// interstimulus interval
      const std::vector<double> rho_;  /// fraction of motor unit types
      const std::vector<double> F_;    /// twitch force
      const std::vector<double> T_;    /// twitch contraction time

      /// stretch dependent activation contribution
      const double lambdaMin_;  /// minimal active fiber stretch
      const double
          lambdaOpt_;  /// optimal active fiber stretch related active nominal stress maximimum

      /// velocity dependent activation contribution
      const double dotLambdaMMin_;  /// minimal stretch rate
      const double ke_;  /// dimensionless constant controlling the curvature of the velocity
                         /// dependent activation function in the eccentric case
      const double kc_;  /// dimensionless constant controlling the curvature of the velocity
                         /// dependent activation function in the concentric case

      /// prescribed activation in corresponding time intervals
      const int actTimesNum_;               /// number of time boundaries
      const std::vector<double> actTimes_;  /// time boundaries between intervals
      const int actIntervalsNum_;           /// number of time intervals
      const std::vector<double>
          actValues_;  /// scaling factor in intervals (1=full activation, 0=no activation)

      /// density
      const double density_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };  // class Muscle_Weickenmeier
  }     // namespace PAR


  class Muscle_WeickenmeierType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "Muscle_WeickenmeierType"; }

    static Muscle_WeickenmeierType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Muscle_WeickenmeierType instance_;
  };  // class Muscle_WeickenmeierType


  class Muscle_Weickenmeier : public So3Material
  {
   public:
    /// Construct empty material object
    Muscle_Weickenmeier();

    /// Construct the material object given material parameters
    explicit Muscle_Weickenmeier(MAT::PAR::Muscle_Weickenmeier* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return Muscle_WeickenmeierType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void Pack(DRT::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// Check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// Material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_muscle_weickenmeier;
    }

    /// Return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Muscle_Weickenmeier(*this));
    }

    /// Material mass density
    double Density() const override { return params_->density_; }

    /// Hyperelastic stress response plus elasticity tensor
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,                ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,                  ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                    ///< Constitutive matrix
        const int gp,                                  ///< Gauss point
        const int eleGID) override;                    ///< Element ID


    /// Setup
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * \brief Post setup routine called before the first Evaluate call to read anisotropy from
     * parameter list
     *
     * \param params Container for additional information
     * \param eleGID Global element id
     */
    void PostSetup(Teuchos::ParameterList& params, int eleGID) override;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
    void EvaluateActiveNominalStress(
        Teuchos::ParameterList& params,  ///< Container for additional information
        const double lambdaM,            ///< Fiber stretch
        double& Pa,                      ///< Active nominal stress
        double& derivPa);  ///< Derivative of active nominal stress w.r.t. the fiber stretch

    /// Evaluate activation level omegaa and its derivative w.r.t. the right Cauchy Green tensor
    void EvaluateActivationLevel(
        Teuchos::ParameterList& params,  ///< Container for additional information
        const double lambdaM,            ///< Fiber stretch
        LINALG::Matrix<3, 3>& M,  ///< Fiber structural tensor (dyadic product of fiber directions)
        double Pa,                ///< Active nominal stress
        double derivPa,           ///< Derivative of active nominal stress w.r.t. the fiber stretch
        double& omegaa,           ///< Activation level
        LINALG::Matrix<3, 3>& domegaadC);  ///< Derivative of the activation level w.r.t. the right
                                           ///< Cauchy Green tensor

    /// Evaluate Lambert W function with Halley's method
    void EvaluateLambert(double xi,  ///< Argument of Lambert W function W(xi)
        double& W0,                  ///< Solution of principal branch of Lambert W function
        double tol,                  ///< Tolerance for Halley's approximation
        double maxiter);             ///< Maximal number of iterations

    /// Add scaled derivative of invC*L*invC to Cmat
    void AddtoCmatDerivInvCLInvCProduct(LINALG::Matrix<6, 6>& cmat,  ///< Constitutive matrix
        const LINALG::Matrix<3, 3>& invC,       ///< Inverse of right Cauchy Green tensor
        const LINALG::Matrix<3, 3>& invCLinvC,  ///< Product invC*L*invC with modified structual
                                                ///< tensor and inverse of right Cauchy Green tensor
        double scalar);                         ///< scaling factor

    /// Weickenmeier material parameters
    MAT::PAR::Muscle_Weickenmeier* params_;

    /// Total simulation time
    double t_tot_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;
  };  // class Muscle_Weickenmeier

}  // namespace MAT

#endif
