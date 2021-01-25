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
      //! @name passive material parameters
      const double alpha_;   ///< material parameter, >0
      const double beta_;    ///< material parameter, >0
      const double gamma_;   ///< material parameter, >0
      const double kappa_;   ///< material parameter for coupled volumetric contribution
      const double omega0_;  ///< weighting factor for isotropic tissue constituents, governs ratio
                             ///< between muscle matrix material (omega0) and muscle fibers (omegap)
                             ///< with omega0 + omegap = 1
                             //! @}

      //! @name active microstructural parameters
      //! @name stimulation frequency dependent activation contribution
      const double
          Na_;  ///< number of active motor units (MU) per undeformed muscle cross-sectional area
      const int muTypesNum_;  ///< number of motor unit (MU) types
                              ///< vectors indicating corresponding parameters of motor unit types
                              ///< e.g. slow type I fibres (index 0), fast resistant  type IIa
                              ///< fibres (index 1), fast fatigue type IIb fibres (index 2)
      const std::vector<double> I_;    ///< interstimulus interval
      const std::vector<double> rho_;  ///< fraction of motor unit types
      const std::vector<double> F_;    ///< twitch force
      const std::vector<double> T_;    ///< twitch contraction time
      //! @}

      //! @name stretch dependent activation contribution
      const double lambdaMin_;  ///< minimal active fiber stretch
      const double
          lambdaOpt_;  ///< optimal active fiber stretch related active nominal stress maximimum
      //! @}

      //! @name velocity dependent activation contribution
      //! @{
      const double dotLambdaMMin_;  ///< minimal stretch rate
      const double ke_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the eccentric case
      const double kc_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the concentric case
      //! @}
      //! @}

      //! @name prescribed activation in corresponding time intervals
      //! @{
      const int actTimesNum_;                ///< number of time boundaries
      const std::vector<double> actTimes_;   ///< time boundaries between intervals
      const int actIntervalsNum_;            ///< number of time intervals
      const std::vector<double> actValues_;  ///< scaling factor in intervals
                                             ///< (1=full activation, 0=no activation)
      //! @}

      /// density
      const double density_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };  // end class Muscle_Weickenmeier
  }     // end namespace PAR


  class Muscle_WeickenmeierType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "Muscle_WeickenmeierType"; }

    static Muscle_WeickenmeierType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Muscle_WeickenmeierType instance_;
  };  // end class Muscle_WeickenmeierType


  class Muscle_Weickenmeier : public So3Material
  {
   public:
    // Construct empty material object
    Muscle_Weickenmeier();

    // Construct the material object given material parameters
    explicit Muscle_Weickenmeier(MAT::PAR::Muscle_Weickenmeier* params);

    // Return unique ParObject id
    int UniqueParObjectId() const override
    {
      return Muscle_WeickenmeierType::Instance().UniqueParObjectId();
    }

    // Pack this class so it can be communicated
    void Pack(DRT::PackBuffer& data) const override;

    // Unpack data from a char vector into this class
    void Unpack(const std::vector<char>& data) override;

    // Check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    // Material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_muscle_weickenmeier;
    }

    // Return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Muscle_Weickenmeier(*this));
    }

    // Material mass density
    double Density() const override { return params_->density_; }

    // Hyperelastic stress response plus elasticity tensor
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  // Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,          // Green-Lagrange strain
        Teuchos::ParameterList& params,                // Container for additional information
        LINALG::Matrix<6, 1>* stress,                  // 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                    // Constitutive matrix
        const int gp,                                  // Gauss point
        const int eleGID) override;                    // Element ID


    // Setup
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    // Post setup routine called before first Evaluate call to read anisotropy from parameter list
    void PostSetup(Teuchos::ParameterList& params, int eleGID) override;

    // Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /*!
     *  \brief Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
     *
     *  \param[in]     params Container for additional information
     *  \param[in]     lambdaM Fiber stretch
     *  \param[in,out] Pa Active nominal stress
     *  \param[in,out] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     */
    /// Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
    void EvaluateActiveNominalStress(
        Teuchos::ParameterList& params, const double lambdaM, double& Pa, double& derivPa);

    /*!
     *  \brief Evaluate activation level omegaa and its derivative w.r.t. the right Cauchy Green
     * tensor
     *
     *  \param[in]     params Container for additional information
     *  \param[in]     lambdaM Fiber stretch
     *  \param[in]     M Fiber structural tensor (dyadic product of fiber directions)
     *  \param[in]     Pa Active nominal stress
     *  \param[in]     derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     *  \param[in,out] omegaa Activation level
     *  \param[in,out] domegaadC Derivative of the activation level w.r.t. the right Cauchy
     *                           Green tensor
     */
    void EvaluateActivationLevel(Teuchos::ParameterList& params, const double lambdaM,
        LINALG::Matrix<3, 3>& M, double Pa, double derivPa, double& omegaa,
        LINALG::Matrix<3, 3>& domegaadC);

    /*!
     *  \brief Evaluate Lambert W function with Halley's method
     *
     *  Solution of Lambert W function is functional inverse of xi = W_0*exp(W_0)
     *  Computation here restricted to principal branch W_0
     *  Use of Halley's method according to:
     *  https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-funsolutiction/
     *
     *  \param[in]     xi Argument of Lambert W function W(xi)
     *  \param[in,out] WO Solution of principal branch of Lambert W function
     *  \param[in]     tol Tolerance for Halley's approximation
     *  \param[in]     maxiter Maximal number of iterations
     */
    void EvaluateLambert(double xi, double& W0, double tol, double maxiter);

    /*!
     *  \brief Add scaled derivative of invC*L*invC to Cmat
     *
     *  Add the following contribution to the constitutive tensor cmat(6,6) based on
     *  - the inverse of the right Cauchy-Green vector invC
     *  - the term invC*L*invC with the structural tensor L
     *
     *              \partial tensor(C)^-1 tensor(L) tensor(C)^-1
     *  scalar *    --------------------------------------------
     *              \partial tensor(C)
     *  wherein the derivative frac{\partial tensor(C)^-1 tensor(L) tensor(C)^-1}{\partial
     *  tensor(C)} is computed to:  - 1/2 * ( Cinv_{ik} Cinv_{jm} L_{mn} Cinv_{nl} + Cinv_{il}
     *  Cinv_{jm} L_{mn} Cinv_{nk} + Cinv_{jk} Cinv_{im} L_{mn} Cinv_{nl} + Cinv_{jl} Cinv_{im}
     *  L_{mn} Cinv_{nk})
     *
     *  \param[in,out] cmat Constitutive matrix
     *  \param[in]     invC Inverse of right Cauchy Green tensor in matrix notation
     *  \param[in]     invCLinvC Product invC*L*invC with modified structual tensor and inverse of
     *                           right Cauchy Green tensor
     *  \param[in]     scalar scaling factor
     */
    void AddtoCmatDerivInvCLInvCProduct(LINALG::Matrix<6, 6>& cmat,
        const LINALG::Matrix<3, 3>& invC, const LINALG::Matrix<3, 3>& invCLinvC, double scalar);

    /// Weickenmeier material parameters
    MAT::PAR::Muscle_Weickenmeier* params_;

    /// Total simulation time
    double t_tot_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;
  };  // end class Muscle_Weickenmeier

}  // end namespace MAT

#endif
