/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Weickenmeier active skeletal muscle material (generalized active strain
approach)


\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MAT_MUSCLE_WEICKENMEIER_H
#define BACI_MAT_MUSCLE_WEICKENMEIER_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_anisotropy.H"
#include "baci_mat_anisotropy_extension_default.H"
#include "baci_mat_anisotropy_extension_provider.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    class Muscle_Weickenmeier : public Parameter
    {
     public:
      /// constructor
      Muscle_Weickenmeier(Teuchos::RCP<MAT::PAR::Material> matdata);

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


  class Muscle_WeickenmeierType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "Muscle_WeickenmeierType"; }

    static Muscle_WeickenmeierType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Muscle_WeickenmeierType instance_;
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
  class Muscle_Weickenmeier : public So3Material
  {
   public:
    // Constructor for empty material object
    Muscle_Weickenmeier();

    // Constructor for the material given the material parameters
    explicit Muscle_Weickenmeier(MAT::PAR::Muscle_Weickenmeier* params);

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Muscle_Weickenmeier(*this));
    }

    MAT::PAR::Parameter* Parameter() const override { return params_; }

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_muscle_weickenmeier;
    };

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    double Density() const override { return params_->density_; }

    int UniqueParObjectId() const override
    {
      return Muscle_WeickenmeierType::Instance().UniqueParObjectId();
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
    MAT::PAR::Muscle_Weickenmeier* params_{};

    /// Fibre stretch of the previous timestep
    double lambdaMOld_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;
  };  // end class Muscle_Weickenmeier

}  // end namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_MUSCLE_WEICKENMEIER_H