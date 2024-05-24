/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Giantesio active skeletal muscle material (active strain approach)

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MUSCLE_GIANTESIO_HPP
#define FOUR_C_MAT_MUSCLE_GIANTESIO_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_four_tensor.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_mat_anisotropy_extension_provider.hpp"
#include "4C_mat_muscle_utils.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_local_numeric_methods.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    class MuscleGiantesio : public CORE::MAT::PAR::Parameter
    {
     public:
      /// constructor
      MuscleGiantesio(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;

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

      //! @name velocity dependent activation contribution: slight modification of original function
      //! given by Weickenmeier et al.
      //! @{
      const double dotLambdaMMin_;  ///< minimal stretch rate
      const double ke_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the eccentric case
      const double kc_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the concentric case
      const double de_;  ///< dimensionless constant controlling the amplitude of the velocity
                         ///< dependent activation function in the eccentric case
      const double dc_;  ///< dimensionless constant controlling the amplitude of the velocity
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

      const double density_;  ///< density
      //@}

    };  // end class Muscle_Giantesio
  }     // end namespace PAR


  class MuscleGiantesioType : public CORE::COMM::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "Muscle_GiantesioType"; }

    static MuscleGiantesioType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MuscleGiantesioType instance_;
  };
  /*!
   * @brief Giantesio muscle material
   *
   * This constituent represents an active hyperelastic muscle material using the active strain
   * approach as described by Giantesio et al.
   *
   * Reference for the material model: G. Giantesio, A. Musesti, 'Strain-dependent internal
   * parameters in hyperelastic biological materials', International Journal of Non-Linear
   * Mechanics, vol. 95, pp. 162-167, 2017, doi:10.1016/j.ijnonlinmec.2017.06.012.
   *
   * Note that the Giantesio muscle material is fully incompressible. Here we introduce a volumetric
   * contribution in the same manner as it is done in Weickenmeier et al. (coupled approach).
   *
   * Reference for the volumetric contribution: J. Weickenmeier, M. Itskov, E Mazza and M. Jabareen,
   * 'A physically motivated constitutive model for 3D numerical simulation of skeletal muscles',
   * International journal for numerical methods in biomedical engineering, vol. 30, no. 5, pp.
   * 545-562, 2014, doi: 10.1002/cnm.2618.
   */
  class MuscleGiantesio : public MAT::So3Material
  {
   public:
    // Constructor for empty material object
    MuscleGiantesio();

    // Constructor for the material given the material parameters
    explicit MuscleGiantesio(MAT::PAR::MuscleGiantesio* params);

    [[nodiscard]] Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new MuscleGiantesio(*this));
    }

    [[nodiscard]] CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

    [[nodiscard]] CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_muscle_giantesio;
    };

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (kinem != INPAR::STR::KinemType::linear && kinem != INPAR::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    [[nodiscard]] double Density() const override { return params_->density_; }

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return MuscleGiantesioType::Instance().UniqueParObjectId();
    }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    bool UsesExtendedUpdate() override { return true; };

    void Update(CORE::LINALG::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params, int const eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

   private:
    /*!
     * @brief Evaluate activation level omegaa and its first and second derivative w.r.t. the fiber
     * stretch
     *
     * The activation level is obtained as the solution of the nonlinear equation given in Giantesio
     * et al. The first and second derivative w.r.t. the fiber stretch is obtained by a central
     * difference approximation.
     *
     * @param[in] lambdaM Fiber stretch
     * @param[in] dotLambdaM Contraction velocity
     * @param[in] currentTime Current time t_n
     * @param[out] omegaaAndDerivs Activation level and derivatives
     */
    CORE::UTILS::ValuesFunctAndFunctDerivs evaluate_activation_level_and_derivatives(
        const double& lambdaM, const double& dotLambdaM, const double& currentTime);

    /*!
     * @brief Solve the nonlinear equation f(omegaa) given in Giantesio et al. for the activation
     * level omegaa.
     *
     * In Giantesio et al. the nonlinear equation is given in the form rhs(omegaa) = lhs(omegaa). To
     * solve for omegaa we solve the equation f(omegaa) = lhs - rhs = 0. This includes the following
     * steps 1) Compute the right hand side 2) Setup the activation level equation as a function lhs
     * - rhs 3) Compute a starting guess for the newton solver using the bisection method 4) Solve
     * the equation using the Newton Raphson method
     *
     * @param[in] lambdaM Fiber stretch
     * @param[in] dotLambdaM Contraction velocity
     * @param[in] currentTime Current time t_n
     * @return omegaa Activation level
     */
    double solve_activation_level_equation(
        const double& lambdaM, const double& dotLambdaM, const double& currentTime);

    /*!
     * @brief Evaluate the activation level equation as given in Giantesio et al. for a given right
     * hand side.
     *
     * In Giantesio et al. the equation is given in the form rhs = lhs. Here we compute the function
     * f(omegaa) = lhs - rhs. The respective derivate w.r.t. the activation level omegaa is
     * df/domegaa = dlhs/domegaa.
     *
     * @param[in] omegaa Activation level
     * @param[in] lambdaM Fiber stretch
     * @param[in] rhs Right hand side of equation
     * @return equationValueAndDeriv Value of the activation level equation and its derivative for
     * the given omegaa, lambdaM and rhs
     */
    std::tuple<double, double> evaluate_activation_level_equation_and_deriv(
        double omegaa, const double& lambdaM, const double& rhs);

    /*!
     * @brief Evaluate the right hand side of the activation level equation as given in Giantesio et
     * al. This includes a time-, stretch- and velocity-dependence.
     *
     * @param[in] lambdaM Fiber stretch
     * @param[in] dotLambdaM Contraction velocity
     * @param[in] currentTime Current time t
     * @return rhs Right hand side of the activation level equation
     */
    double evaluate_rhs_activation_level_equation(
        const double& lambdaM, const double& dotLambdaM, const double& currentTime);

    /*!
     * @brief Evaluate the integral of the active nominal stress w.r.t. the fiber stretch as given
     * in Giantesio et al. This includes a time-, stretch- and velocity-dependence.
     *
     * @param[in] lambdaM Fiber stretch
     * @param[in] dotLambdaM Contraction velocity
     * @param[in] currentTime Current time t
     * @return intPa Integral of the active nominal stress w.r.t. the fiber stretch
     */
    double evaluate_active_nominal_stress_integral(
        const double& lambdaM, const double& dotLambdaM, const double& currentTime);

    /*!
     * @brief Returns the active deformation gradient Fa
     *
     * @param[in] omegaa Activation level
     * @param[in] M Structural tensor of fiber directions
     * @return Fa Active deformation gradient
     */
    CORE::LINALG::Matrix<3, 3> act_def_grad(
        const double omegaa, const CORE::LINALG::Matrix<3, 3>& M);

    /*!
     * @brief Returns the first derivative of the active deformation gradient Fa w.r.t. omegaa
     *
     * @param[in] omegaa Activation level
     * @param[in] M Structural tensor of fiber directions
     * @return dFadomegaa First derivative of the active deformation gradient
     */
    CORE::LINALG::Matrix<3, 3> d_act_def_grad_d_act_level(
        const double omegaa, const CORE::LINALG::Matrix<3, 3>& M);

    /*!
     * @brief Returns the second derivative of the active deformation gradient Fa w.r.t. omegaa
     *
     * @param[in] omegaa Activation level
     * @param[in] M Structural tensor of fiber directions
     * @return ddFaddomegaa Second derivative of the active deformation gradient
     */
    CORE::LINALG::Matrix<3, 3> dd_act_def_grad_dd_act_level(
        const double omegaa, const CORE::LINALG::Matrix<3, 3>& M);

    /*!
     * @brief Returns the inverse of the active deformation gradient Fa
     *
     * @param[in] Fa Active deformation gradient
     * @return invFa Inverse of the active deformation gradient
     */
    CORE::LINALG::Matrix<3, 3> inv_act_def_grad(const CORE::LINALG::Matrix<3, 3>& Fa);

    /*!
     * @brief Returns the first derivative of the inverse of the active deformation gradient Fa
     * w.r.t. omegaa
     *
     * @param[in] Fa Active deformation gradient
     * @param[in] dFadomegaa First derivative of the active deformation gradient
     * @return dinvFadomegaa First derivative of the inverse of the active deformation gradient
     */
    CORE::LINALG::Matrix<3, 3> d_inv_act_def_grad_d_act_level(
        const CORE::LINALG::Matrix<3, 3>& Fa, const CORE::LINALG::Matrix<3, 3>& dFadomegaa);

    /*!
     * @brief Check if material is active at current time and fiber stretch lambdaM
     *
     * @param[in] currentTime Active deformation gradient
     * @return isActive true, if active; false, if passive
     */
    bool is_active(const double& currentTime);

    /// Giantesio material parameters
    MAT::PAR::MuscleGiantesio* params_{};

    /// Fibre stretch of the previous timestep
    double lambda_m_old_;

    /// Activation level of the previous timestep
    double omegaa_old_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropy_extension_;
  };  // end class Muscle_Giantesio

}  // end namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
