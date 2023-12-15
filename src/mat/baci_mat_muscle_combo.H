/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Combo active skeletal muscle material (modified and corrected generalized
active strain approach) with variable time-dependent activation

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MAT_MUSCLE_COMBO_H
#define BACI_MAT_MUSCLE_COMBO_H

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_anisotropy.H"
#include "baci_mat_anisotropy_extension_default.H"
#include "baci_mat_anisotropy_extension_provider.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"
#include "baci_utils_function.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    class Muscle_Combo : public Parameter
    {
     public:
      /// constructor
      Muscle_Combo(Teuchos::RCP<MAT::PAR::Material> matdata);

      Teuchos::RCP<MAT::Material> CreateMaterial() override;

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

      //! @name function-id in input file for user-prescribed time-dependent activation
      const int actFunctId_;

      const double density_;  ///< density
      //@}

    };  // end class Muscle_Combo
  }     // end namespace PAR


  class Muscle_ComboType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "Muscle_ComboType"; }

    static Muscle_ComboType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Muscle_ComboType instance_;
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
  class Muscle_Combo : public So3Material
  {
   public:
    // Constructor for empty material object
    Muscle_Combo();

    // Constructor for the material given the material parameters
    explicit Muscle_Combo(MAT::PAR::Muscle_Combo* params);

    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Muscle_Combo(*this)); }

    MAT::PAR::Parameter* Parameter() const override { return params_; }

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_muscle_combo; };

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    double Density() const override { return params_->density_; }

    int UniqueParObjectId() const override
    {
      return Muscle_ComboType::Instance().UniqueParObjectId();
    }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    bool UsesExtendedUpdate() override { return true; };

    void Update(CORE::LINALG::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params, int const eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

   private:
    /*!
     *  \brief Evaluate active nominal stress Pa, its integral and its derivative w.r.t. the fiber
     * stretch
     *
     *  \param[in] params Container for additional information
     *  \param[in] lambdaM Fiber stretch
     *  \param[in] intPa Integral of the active nominal stress from lambdaMin to lambdaM
     *  \param[out] Pa Active nominal stress
     *  \param[out] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     */
    void EvaluateActiveNominalStress(Teuchos::ParameterList& params, const double lambdaM,
        double& intPa, double& Pa, double& derivPa);

    /*!
     *  \brief Evaluate activation level omegaa and its first and second derivatives w.r.t. the
     * fiber stretch
     *
     *  \param[in] lambdaM Fiber stretch
     *  \param[in] intPa Integral of the active nominal stress from lambdaMin to lambdaM
     *  \param[in] Pa Active nominal stress
     *  \param[in] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     *  \param[out] omegaa Activation level
     *  \param[out] derivOmegaa Derivative of the activation level w.r.t. the fiber stretch
     *  \param[out] derivDerivOmegaa Second derivative of the activation level w.r.t. the fiber
     * stretch
     */
    void EvaluateActivationLevel(const double lambdaM, const double intPa, const double Pa,
        const double derivPa, double& omegaa, double& derivOmegaa, double& derivDerivOmegaa);

    /// Combo material parameters
    MAT::PAR::Muscle_Combo* params_{};

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;

    /// Activation function
    const CORE::UTILS::FunctionOfSpaceTime* activationFunction_;
  };  // end class Muscle_Combo

}  // end namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_MUSCLE_COMBO_H