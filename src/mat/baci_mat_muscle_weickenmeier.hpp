/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Weickenmeier active skeletal muscle material (generalized active strain
approach)


\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MUSCLE_WEICKENMEIER_HPP
#define FOUR_C_MAT_MUSCLE_WEICKENMEIER_HPP

#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_anisotropy.hpp"
#include "baci_mat_anisotropy_extension_default.hpp"
#include "baci_mat_anisotropy_extension_provider.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    class MuscleWeickenmeier : public Parameter
    {
     public:
      /// constructor
      MuscleWeickenmeier(Teuchos::RCP<MAT::PAR::Material> matdata);

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
      //! given by Weickenmeier et al..
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

    };  // end class Muscle_Weickenmeier
  }     // end namespace PAR


  class MuscleWeickenmeierType : public CORE::COMM::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "Muscle_WeickenmeierType"; }

    static MuscleWeickenmeierType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MuscleWeickenmeierType instance_;
  };


  /*!
   * \brief Weickenmeier muscle material
   *
   * This constituent represents an active hyperelastic muscle material using the generalized active
   * strain approach as described by Weickenmeier et al.
   *
   * Reference for the material model: J. Weickenmeier, M. Itskov, E Mazza and M. Jabareen, 'A
   * physically motivated constitutive model for 3D numerical simulation of skeletal muscles',
   * International journal for numerical methods in biomedical engineering, vol. 30, no. 5, pp.
   * 545-562, 2014, doi: 10.1002/cnm.2618.
   *
   * Amongst others, the Weickenmeier material accounts for the velocity dependence of the active
   * force production. However, the given function produces unsuitable values - in the sense that it
   * fails to reproduce experimetally observerd force-velocity curves. Thus, a slight modification
   * of the velocity dependence in the style of Boel et al.. is introduced here. Boel et al.. scale
   * and shift the eccentric branch of the velocity dependency using the parameter d. To be able to
   * turn the velocity dependency on and off, this concept is adopted for the concentric branch as
   * well. Choosing de=dc=0 leads to a neglection of the velocity dependence. Choosing dc=1 and
   * de=d-1 reproduces the function presented by Boel et al..
   *
   * Reference: M. Boel and S. Reese, 'Micromechanical modelling of skeletal muscles based on the
   * finite element method', Computer Methods in Biomechanics and Biomedical Engineering, vol. 11,
   * no. 5, pp. 489-504, 2008, doi: 10.1080/10255840701771750.
   */
  class MuscleWeickenmeier : public So3Material
  {
   public:
    // Constructor for empty material object
    MuscleWeickenmeier();

    // Constructor for the material given the material parameters
    explicit MuscleWeickenmeier(MAT::PAR::MuscleWeickenmeier* params);

    [[nodiscard]] Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new MuscleWeickenmeier(*this));
    }

    [[nodiscard]] MAT::PAR::Parameter* Parameter() const override { return params_; }

    [[nodiscard]] INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_muscle_weickenmeier;
    };

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (kinem != INPAR::STR::KinemType::linear && kinem != INPAR::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    [[nodiscard]] double Density() const override { return params_->density_; }

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return MuscleWeickenmeierType::Instance().UniqueParObjectId();
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
     *  \brief Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
     *
     *  \param[in] params Container for additional information
     *  \param[in] lambdaM Fiber stretch
     *  \param[out] Pa Active nominal stress
     *  \param[out] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     */
    void EvaluateActiveNominalStress(
        Teuchos::ParameterList& params, const double lambdaM, double& Pa, double& derivPa);

    /*!
     *  \brief Evaluate activation level omegaa and its derivative w.r.t. the fiber stretch
     *
     *  \param[in] params Container for additional information
     *  \param[in] lambdaM Fiber stretch
     *  \param[in] Pa Active nominal stress
     *  \param[in] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     *  \param[out] omegaa Activation level
     *  \param[out] derivOmegaa Derivative of the activation level w.r.t. the fiber stretch
     */
    void EvaluateActivationLevel(const double lambdaM, const double Pa, const double derivPa,
        double& omegaa, double& derivOmegaa);

    /// Weickenmeier material parameters
    MAT::PAR::MuscleWeickenmeier* params_{};

    /// Fibre stretch of the previous timestep
    double lambda_m_old_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropy_extension_;
  };  // end class Muscle_Weickenmeier

}  // end namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif