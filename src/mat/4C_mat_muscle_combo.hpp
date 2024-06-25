/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Combo active skeletal muscle material (modified and corrected generalized
active strain approach) with variable time-dependent activation

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MUSCLE_COMBO_HPP
#define FOUR_C_MAT_MUSCLE_COMBO_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_material.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_RCP.hpp>

#include <unordered_map>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    class MuscleCombo : public Core::Mat::PAR::Parameter
    {
     public:
      /// constructor
      MuscleCombo(const Core::Mat::PAR::Parameter::Data& matdata);

      Teuchos::RCP<Core::Mat::Material> create_material() override;

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
      const double Popt_;  ///< optimal (maximal) active tetanised stress
      //! @}

      //! @name stretch dependent activation contribution
      const double lambdaMin_;  ///< minimal active fiber stretch
      const double
          lambdaOpt_;  ///< optimal active fiber stretch related active nominal stress maximimum
      //! @}

      //! @name time-/space-dependent activation

      //! type of activation prescription
      const Inpar::Mat::ActivationType activationType_;

      /*!
       * @brief type-dependent parameters for activation
       *
       * Depending on the type of activation prescription this is one of the options below:
       * - Id of the function in the input file specifying an analytical function
       * - Map retrieved from the pattern file path in the input file specifying a discrete values.
       *   The integer key refers to the elememt ids, the vector bundles time-activation pairs.
       */
      using ActivationParameterVariant = std::variant<std::monostate, const int,
          const std::unordered_map<int, std::vector<std::pair<double, double>>>>;
      ActivationParameterVariant activationParams_;
      //! @}

      const double density_;  ///< density
      //@}

    };  // end class MuscleCombo
  }     // end namespace PAR


  class MuscleComboType : public Core::Communication::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "Muscle_ComboType"; }

    static MuscleComboType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MuscleComboType instance_;
  };


  /*!
   * \brief Combo muscle material
   *
   * This constituent represents an active hyperelastic muscle material using a generalized active
   * strain approach. Stress and material tangent are consistently derived from the strain-energy
   * function.
   *
   * The general material formulation is equal to Weickenmeier et al. [1] with the following
   * modifications:
   * 1. The derivative of omegaa w.r.t. C is included as described for the active stress
   * approach in Giantesio et al. [2]. This leads to an additional term in the stress and material
   * tangent computation and an updated equation for the activation level omegaa.
   * 2. The twitch superposition is neglected and the time and space depenent optimal nominal stress
   * is computed through a user-prescribed function.
   * 3. A velocity dependence is not considered.
   *
   * References:
   * [1] J. Weickenmeier, M. Itskov, E Mazza and M. Jabareen, 'A
   * physically motivated constitutive model for 3D numerical simulation of skeletal muscles',
   * International journal for numerical methods in biomedical engineering, vol. 30, no. 5, pp.
   * 545-562, 2014, doi: 10.1002/cnm.2618.
   * [2] G. Giantesio, A. Musesti, 'Strain-dependent internal
   * parameters in hyperelastic biological materials', International Journal of Non-Linear
   * Mechanics, vol. 95, pp. 162-167, 2017, doi:10.1016/j.ijnonlinmec.2017.06.012.
   */
  class MuscleCombo : public So3Material
  {
   public:
    // Constructor for empty material object
    MuscleCombo();

    // Constructor for the material given the material parameters
    explicit MuscleCombo(Mat::PAR::MuscleCombo* params);

    [[nodiscard]] Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new MuscleCombo(*this));
    }

    [[nodiscard]] Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    [[nodiscard]] Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_muscle_combo;
    };

    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (kinem != Inpar::STR::KinemType::linear && kinem != Inpar::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    [[nodiscard]] double Density() const override { return params_->density_; }

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return MuscleComboType::Instance().UniqueParObjectId();
    }

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(const std::vector<char>& data) override;

    void setup(int numgp, Input::LineDefinition* linedef) override;

    bool UsesExtendedUpdate() override { return true; };

    void update(Core::LinAlg::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params, int const eleGID) override;

    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    using ActivationEvaluatorVariant =
        std::variant<std::monostate, const Core::UTILS::FunctionOfSpaceTime*,
            const std::unordered_map<int, std::vector<std::pair<double, double>>>*>;

   private:
    /*!
     * \brief Evaluate active nominal stress Pa, its integral and its derivative w.r.t. the fiber
     * stretch
     *
     * \param[in] params Container for additional information
     * \param[in] eleGID Global element id used for discretely prescribed activation
     * \param[in] lambdaM Fiber stretch
     * \param[in] intPa Integral of the active nominal stress from lambdaMin to lambdaM
     * \param[out] Pa Active nominal stress
     * \param[out] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     */
    void evaluate_active_nominal_stress(Teuchos::ParameterList& params, const int eleGID,
        const double lambdaM, double& intPa, double& Pa, double& derivPa);

    /*!
     * \brief Evaluate activation level omegaa and its first and second derivatives w.r.t. the
     * fiber stretch
     *
     * \param[in] lambdaM Fiber stretch
     * \param[in] intPa Integral of the active nominal stress from lambdaMin to lambdaM
     * \param[in] Pa Active nominal stress
     * \param[in] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     * \param[out] omegaa Activation level
     * \param[out] derivOmegaa Derivative of the activation level w.r.t. the fiber stretch
     * \param[out] derivDerivOmegaa Second derivative of the activation level w.r.t. the fiber
     * stretch
     */
    void evaluate_activation_level(const double lambdaM, const double intPa, const double Pa,
        const double derivPa, double& omegaa, double& derivOmegaa, double& derivDerivOmegaa);

    /// Combo material parameters
    Mat::PAR::MuscleCombo* params_{};

    /// Holder for anisotropic behavior
    Mat::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    Mat::DefaultAnisotropyExtension<1> anisotropy_extension_;

    /// Activation evaluator, either analytical symbolic function of space and time or discrete
    /// activation map
    ActivationEvaluatorVariant activation_evaluator_;
  };  // end class Muscle_Combo

}  // end namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif