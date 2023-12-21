/*--------------------------------------------------------------------------*/
/*! \file
\brief An elastic visco-plastic finite strain material law without yield surface.

\level 2
*--------------------------------------------------------------------------*/

#ifndef BACI_MAT_VISCOPLASTIC_NO_YIELD_SURFACE_H
#define BACI_MAT_VISCOPLASTIC_NO_YIELD_SURFACE_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for the ViscoPlasticNoYieldSurface material
    class ViscoPlasticNoYieldSurface : public Parameter
    {
     public:
      explicit ViscoPlasticNoYieldSurface(Teuchos::RCP<MAT::PAR::Material> matdata);
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! @name return methods of material parameters
      //! @{

      //! return Density  \f$ \rho \f$
      double Density() const { return density_; };
      //! return Poisson's ratio \f$ \nu \f$
      double Nue() const { return nue_; };
      //! return Young's modulus  \f$ E \f$
      double Young() const { return young_; };
      //! return temperature \f$ T \f$
      double Temperature() const { return temperature_; };
      //! return plastic shear strain rate pre-exp-factor \f$ A \f$
      double PreExpFac() const { return pre_exp_fac_; };
      //! return activation energy \f$ Q \f$
      double ActivationEnergy() const { return activation_energy_; };
      //! return gas constant \f$ R \f$
      double GasConstant() const { return gas_constant_; };
      //! return strain-rate-sensitivity \f$ m \f$
      double StrainRateSensitivity() const { return strain_rate_sensitivity_; };
      //! return flow resistance pre-factor \f$ H_0 \f$
      double FlowResPreFac() const { return flow_res_pre_fac_; };
      //! return initial flow resistance \f$ S^0 \f$
      double InitFlowRes() const { return init_flow_res_; };
      //! return flow resistance exponent \f$ a \f$
      double FlowResExp() const { return flow_res_exp_; };
      //! return flow resistance saturation factor \f$ S_* \f$
      double FlowResSatFac() const { return flow_res_sat_fac_; };
      //! return flow resistance saturation exponent \f$ b \f$
      double FlowResSatExp() const { return flow_res_sat_exp_; };

      //! @}

     private:
      //! @name material parameters
      //! @{

      //! Density \f$ \rho \f$
      const double density_;
      //! Possion's ratio \f$ \nu \f$
      const double nue_;
      //! Young's modulus \f$ E \f$
      const double young_;
      //! temperature \f$ T \f$
      const double temperature_;
      //! plastic shear strain rate pre-exp-factor \f$ A \f$
      const double pre_exp_fac_;
      //! activation energy \f$ Q \f$
      const double activation_energy_;
      //! gas constant \f$ R \f$
      const double gas_constant_;
      //! strain-rate-sensitivity \f$ m \f$
      const double strain_rate_sensitivity_;
      //! initial flow resistance \f$ S^0 \f$
      const double init_flow_res_;
      //! flow resistance exponent \f$ a \f$
      const double flow_res_exp_;
      //! flow resistance pre-factor \f$ H_0 \f$
      const double flow_res_pre_fac_;
      //! flow resistance saturation factor \f$ S_* \f$
      const double flow_res_sat_fac_;
      //! flow resistance saturation exponent \f$ b \f$
      const double flow_res_sat_exp_;

      //! @}
    };
  }  // namespace PAR

  class ViscoPlasticNoYieldSurfaceType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ViscoPlasticNoYieldSurfaceType"; }

    static ViscoPlasticNoYieldSurfaceType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ViscoPlasticNoYieldSurfaceType instance_;
  };

  /*----------------------------------------------------------------------*/
  // forward declaration
  struct PreCalculatedTerms;

  /*!
   * \brief Implementation of a visco-plastic finite strain material with no yield surface.
   *
   * This class implements the visco-plastic finite strain material with no yield surface based on
   * the following papers / student work:
   *   -# L. Anand et al., An Elastic-Viscoplastic Model for Lithium", Journal of The
   *      Electrochemical Society, 2019.
   *   -# G. Weber et al., Finite deformation constitutive equations and a time integration
   *      procedure for isotropic, hyperelastic-viscoplastic solids.", Computer Methods in Applied
   *      Mechanics and Engineering, 1990.
   *   -# Details on the implementation in master's thesis of Philipp Bofinger, Supervisor:
   *      Christoph Schmidt.
   *
   * \note The time integration of the internal evolution equations is a backward one-step time
   *       integration algorithm.
   */
  class ViscoPlasticNoYieldSurface : public So3Material
  {
   public:
    //! construct empty material object
    ViscoPlasticNoYieldSurface();

    //! construct the material object given material parameters
    explicit ViscoPlasticNoYieldSurface(MAT::PAR::ViscoPlasticNoYieldSurface* params);

    int UniqueParObjectId() const override
    {
      return ViscoPlasticNoYieldSurfaceType::Instance().UniqueParObjectId();
    }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

   private:
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_vp_no_yield_surface;
    }

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (kinem != INPAR::STR::kinem_nonlinearTotLag)
        dserror("element and material kinematics are not compatible");
    }

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ViscoPlasticNoYieldSurface(*this));
    }

    //! @name Evaluation methods
    //! @{

    /*!
     * @param[in] Me_trial_Vstress  trial stresses in stress-like Voigt notation
     * @param[in] p                 mean normal pressure
     * @return deviatoric trial stresses
     */
    CORE::LINALG::Matrix<3, 3>& CalculateDeviatoricTrialStresses(
        const CORE::LINALG::Matrix<6, 1>& Me_trial_Vstress, const double p) const;

    /*!
     * @brief elastic stiffness tensor in intermediate configuration is calculated
     *
     * @param[in] eigen_vectors  eigenvectors of \f$ F^{*}_\text{e} \f$
     * @param[in] eigen_values   eigenvalues of \f$ F^{*}_\text{e} \f$
     * @return elastic stiffness tensor in intermediate configuration
     */
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& CalculateElasticStiffness(
        const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
        const CORE::LINALG::Matrix<3, 1>& eigen_values) const;

    /*!
     * @brief calculate linearization for local newton loop of internal evolution equations
     *
     * @param[in] equ_tens_stress_np  current equivalent tensile stress
     * @param[in] flow_resistance_np  current flow resistance
     * @param[in] terms               terms that can be precalculated and reused several times
     * @return Linearization matrix
     */
    CORE::LINALG::Matrix<2, 2>& CalculateLinearization(const double equ_tens_stress_np,
        const double flow_resistance_np, const MAT::PreCalculatedTerms& terms);

    /*!
     * @brief calculate the residual for the equations solved within the local newton loop of
     * internal evolution equations
     *
     * @param[in] equ_tens_stress_np     current equivalent tensile stress
     * @param[in] equ_tens_trial_stress  last converged equivalent tensile stress
     * @param[in] flow_resistance_np     current flow resistance
     * @param[in] flow_resistance_n      last converged flow resistance
     * @param[in] terms                  terms that can be precalculated and reused several times
     * @return residual vector
     */
    CORE::LINALG::Matrix<2, 1>& CalculateResidual(const double equ_tens_stress_np,
        const double equ_tens_trial_stress, const double flow_resistance_np,
        const double flow_resistance_n, const PreCalculatedTerms& terms);

    /*!
     * @param[in] defgrd    deformation gradient
     * @param[in] Re_trial  trial elastic rotation tensor
     * @param[in] Me        energy conjugated stress tensor
     * @return second Piola--Kirchhoff stresses
     */
    CORE::LINALG::Matrix<3, 3>& CalculateSecondPiolaKirchhoffStresses(
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 3>& Re_trial,
        const CORE::LINALG::Matrix<3, 3>& Me) const;

    /*!
     * @param[in] Fe_trial        trial elastic deformation gradient \f$ F^{*}_\text{e} \f$
     * @param[out] eigen_values   eigenvalues of \f$ F^{*}_\text{e} \f$
     * @param[out] eigen_vectors  eigenvectors of \f$ F^{*}_\text{e} \f$
     */
    void CalculateTrialElasticDefgradEigenvaluesAndEigenvectors(
        const CORE::LINALG::Matrix<3, 3>& Fe_trial, CORE::LINALG::Matrix<3, 1>& eigen_values,
        CORE::LINALG::Matrix<3, 3>& eigen_vectors) const;

    /*!
     * @param[in] Fe_trial       trial elastic deformation gradient \f$ F^{*}_\text{e} \f$
     * @param[in] eigen_vectors  eigenvectors of \f$ F^{*}_\text{e} \f$
     * @param[in] eigen_values   eigenvalues of \f$ F^{*}_\text{e} \f$
     * @return trial elastic rotation tensor \f$ R^{*}_\text{e} \f$
     */
    CORE::LINALG::Matrix<3, 3>& CalculateTrialElasticRotation(
        const CORE::LINALG::Matrix<3, 3>& Fe_trial, const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
        const CORE::LINALG::Matrix<3, 1>& eigen_values) const;

    /*!
     * @param[in] eigen_vectors  eigenvectors of \f$ F^{*}_\text{e} \f$
     * @param[in] eigen_values   eigenvalues of \f$ F^{*}_\text{e} \f$
     * @return logarithmic elastic strains in strain-like Voigt notation
     */
    CORE::LINALG::Matrix<6, 1>& CalculateLogElasticStrainInStrainLikeVoigtNotation(
        const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
        const CORE::LINALG::Matrix<3, 1>& eigen_values) const;

    /*!
     * @param[in] Me_trial_dev  deviatoric part of trial elastic stresses
     * @return trial elastic equivalent stress [sqrt(3/2 * trace(Me_trial_dev . Me_trial_dev))]
     */
    double CalculateTrialEquivalentStress(const CORE::LINALG::Matrix<3, 3>& Me_trial_dev) const;

    /*!
     * @brief inverse viscous deformation gradient is updated and returned
     *
     * @param[in] last_iFv       last converged inverse viscous deformation gradient
     * @param[in] eigen_vectors  eigenvectors of \f$ F^{*}_\text{e} \f$
     * @param[in] eigen_values   eigenvalues of \f$ F^{*}_\text{e} \f$
     * @param[in] eta            deviatoric elastic deformation share
     * @return inverse viscous deformation gradient
     */
    CORE::LINALG::Matrix<3, 3>& CalculateUpdatedInverseViscousDefgrad(
        const CORE::LINALG::Matrix<3, 3>& last_iFv, const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
        const CORE::LINALG::Matrix<3, 1>& eigen_values, const double eta) const;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    /*!
     * @brief Local Newton-loop to solve for the internal evolution equations
     *
     * @param[in/out] x  solution vector
     * @param[in]     dt time step
     */
    void LocalNewtonLoop(CORE::LINALG::Matrix<2, 1>& x, double dt);

    /*!
     * @brief Calculate and return terms of the formulation that can be reused several times
     *
     * @param[in] equ_tens_stress_np  current equivalent tensile stress
     * @param[in] flow_resistance_np  current flow resistance
     * @param[in] dt                  time step
     * @return pre-calculated terms
     */
    PreCalculatedTerms PreCalculateTerms(
        const double equ_tens_stress_np, const double flow_resistance_np, const double dt);

    //! @}

    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /*!
     * @brief computes isotropic elasticity tensor in matrix notion for 3d
     *
     * @param[out] cmat elasticity tensor in intermediate configuration
     */
    void SetupCmat(CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat);

    void Update() override;

    bool NeedsDefgrd() override { return true; };

    MAT::PAR::Parameter* Parameter() const override { return params_; }

    double Density() const override { return params_->Density(); }

    //! material parameters
    MAT::PAR::ViscoPlasticNoYieldSurface* params_;

    //! inverse plastic deformation gradient for each Gauss point at last converged state
    std::vector<CORE::LINALG::Matrix<3, 3>> last_plastic_defgrd_inverse_;

    //! current inverse plastic deformation gradient for each Gauss point
    std::vector<CORE::LINALG::Matrix<3, 3>> current_plastic_defgrd_inverse_;

    //! flow resistance 'S' for each Gauss point at last converged state
    std::vector<double> last_flowres_isotropic_;

    //! current flow resistance 'S' for each Gauss point
    std::vector<double> current_flowres_isotropic_;
  };
}  // namespace MAT

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // MAT_VISCOPLASTIC_NO_YIELD_SURFACE_H
