// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOPLASTIC_LAWS_HPP
#define FOUR_C_MAT_VISCOPLASTIC_LAWS_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

#include <memory>


FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class ViscoplasticLawReformulatedJohnsonCook
     *
     * Parameter class of ViscoplasticLawReformulatedJohnsonCook
     */
    class ViscoplasticLawReformulatedJohnsonCook : public Core::Mat::PAR::Parameter
    {
     public:
      explicit ViscoplasticLawReformulatedJohnsonCook(
          const Core::Mat::PAR::Parameter::Data& matdata);

      std::shared_ptr<Core::Mat::Material> create_material() override { return nullptr; };

      // getter methods
      //! get strain rate prefactor \f$ \dot{P}_0 \f$
      [[nodiscard]] double strain_rate_pre_fac() const { return strain_rate_prefac_; };
      //! get exponential factor \f$ C \f$
      [[nodiscard]] double strain_rate_exp_fac() const { return strain_rate_exp_fac_; };
      //! get initial yield strength (\f$ A_0 \f$ scaled by the temperature dependence factor)
      [[nodiscard]] double init_yield_strength() const { return init_yield_strength_; };
      //! get prefactor of the isotropic hardening stress (\f$ B_0 \f$ scaled by the temperature
      //! dependence factor)
      [[nodiscard]] double isotrop_harden_prefac() const { return isotrop_harden_prefac_; };
      //! get exponent of the isotropic hardening stress \f$ n \f$
      [[nodiscard]] double isotrop_harden_exp() const { return isotrop_harden_exp_; };

     private:
      //! strain rate prefactor \f$ \dot{P}_0 \f$
      const double strain_rate_prefac_;

      //! exponential factor \f$ C \f$
      const double strain_rate_exp_fac_;

      //! initial yield strength (\f$ A_0 \f$ scaled by the temperature dependence factor)
      double init_yield_strength_;

      //! prefactor of the isotropic hardening stress (\f$ B_0 \f$ scaled by the temperature
      //! dependence factor)
      double isotrop_harden_prefac_;

      //! exponent of the isotropic hardening stress \f$ n \f$
      const double isotrop_harden_exp_;
    };
  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class ViscoplasticLaws
   * \brief Implementation of an interface called by the viscoplasticity finite strain model
   * InelasticDefgradTransvIsotropElastViscoplast to evaluate the viscoplastic flow rule (plastic
   * strain rate in terms of stress and plasticity state) and the hardening law (increase in yield
   * stress / flow resistance).
   *
   * This class implements an interface for viscoplasticity laws (combinations of flow rule and
   * hardening model), used by the (transversely) isotropic
   * material model InelasticDefgradTransvIsotropElastViscoplast to describe viscoplastic material
   * response. Specifically, it is used to evaluate the plastic strain rate for the given local
   * stress, taking into account the current plasticity state of the material (current plastic
   * strain). It also provides the linearization of this quantity with respect to both the stress
   * and the plastic strain, required in the linearization of the material referenced above. Both
   * viscoplasticity laws with and without yield surfaces can be implemented.
   *
   * For further information on the implemented interface, refer to :
   *    -# Master's Thesis : Dragos-Corneliu Ana, Continuum Modeling and Calibration of
   * Viscoplasticity in the Context of the Lithium Anode in Solid State Batteries, Supervisor:
   * Christoph Schmidt, 2024
   */
  class ViscoplasticLaws
  {
   public:
    /// construct viscoplastic laws with specific material params
    explicit ViscoplasticLaws(Core::Mat::PAR::Parameter* params);
    /// construct empty viscoplastic law
    ViscoplasticLaws();

    /// virtual destructor
    virtual ~ViscoplasticLaws() = default;

    /*!
     * @brief create object by input parameter ID
     *
     * @param[in] matnum  material ID
     * @return pointer to material that is defined by material ID
     */
    static std::shared_ptr<ViscoplasticLaws> factory(int matnum);

    /// provide material type
    virtual Core::Materials::MaterialType material_type() const = 0;

    /*!
     * @brief Evaluate the ratio of the equivalent stress to the yield stress (or for flow rules
     * without yield surface: ratio of the equivalent stress to the flow resistance)
     *
     * @param[in] equiv_stress Equivalent stress
     * @param[in] equiv_plastic_strain Equivalent plastic strain
     * @return Stress ratio: equiv_stress/yield_stress(equiv_plastic_strain)
     */
    virtual double evaluate_stress_ratio(
        const double equiv_stress, const double equiv_plastic_strain) = 0;

    /*!
     * @brief Evaluate the equivalent plastic strain rate for a given equivalent stress and a given
     * plastic strain
     *
     * @param[in] equiv_stress Equivalent stress
     * @param[in] equiv_plastic_strain Equivalent plastic strain
     * @param[in] dt Time step size (used solely for overflow checking purposes)
     * @return Equivalent plastic strain rate
     */
    virtual double evaluate_plastic_strain_rate(const double equiv_stress,
        const double equiv_plastic_strain, const double dt, const bool update_hist_var = true) = 0;

    /*!
     * @brief Evaluate the derivatives of the equivalent plastic strain rate for a given equivalent
     * stress and a given plastic strain
     *
     * @param[in] equiv_stress Equivalent stress
     * @param[in] equiv_plastic_strain Equivalent plastic strain
     * @param[in] dt Time step size (used solely for overflow checking purposes)
     * @return Derivatives of the equivalent plastic strain rate w.r.t. the equivalent stress
     *         (element 0 of matrix) and the plastic strain (element 1 of matrix)
     */
    virtual Core::LinAlg::Matrix<2, 1> evaluate_derivatives_of_plastic_strain_rate(
        const double equiv_stress, const double equiv_plastic_strain, const double dt,
        const bool update_hist_var = true) = 0;

    /// Return material parameters
    virtual Core::Mat::PAR::Parameter* parameter() const { return params_; }

    /*!
     * @brief Setup viscoplasticity law for the specific element
     *
     * @param[in] numgp Number of Gauss points
     * @param[in] container Input parameter Container
     */
    virtual void setup(const int numgp, const Core::IO::InputParameterContainer& container) = 0;

    /*!
     * @brief Pre-evaluation, intended to be used for stuff that has to be done only once per
     *        evaluate()
     *
     * @param[in] gp      Current Gauss point
     */
    virtual void pre_evaluate(int gp) = 0;

    /*!
     * @brief Update history variables of the viscoplasticity law for next time step
     */
    virtual void update() = 0;

    /*!
     * @brief Update the history variables for a specific GP after a converged substep (substepping
     * procedure)
     *
     * @param[in] gp      Gauss point
     */
    virtual void update_gp_state(int gp) = 0;

    virtual void pack_viscoplastic_law(Core::Communication::PackBuffer& data) const = 0;

    virtual void unpack_viscoplastic_law(Core::Communication::UnpackBuffer& buffer) = 0;

   private:
    /// material parameters
    Core::Mat::PAR::Parameter* params_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class ViscoplasticLawReformulatedJohnsonCook
   *
   *  Viscoplasticity law associated to the Reformulated Johnson-Cook model proposed in
   *  -# Mareau et al., A thermodynamically consistent formulation of the Johnson-Cook model,
   *     Mechanics of Materials 143, 2020
   */
  class ViscoplasticLawReformulatedJohnsonCook : public ViscoplasticLaws
  {
   public:
    explicit ViscoplasticLawReformulatedJohnsonCook(Core::Mat::PAR::Parameter* params);

    Mat::PAR::ViscoplasticLawReformulatedJohnsonCook* parameter() const override
    {
      return dynamic_cast<Mat::PAR::ViscoplasticLawReformulatedJohnsonCook*>(
          Mat::ViscoplasticLaws::parameter());
    }

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::mvl_reformulated_Johnson_Cook;
    };

    double evaluate_stress_ratio(
        const double equiv_stress, const double equiv_plastic_strain) override;

    double evaluate_plastic_strain_rate(const double equiv_stress,
        const double equiv_plastic_strain, const double dt, const bool update_hist_var) override;

    Core::LinAlg::Matrix<2, 1> evaluate_derivatives_of_plastic_strain_rate(
        const double equiv_stress, const double equiv_plastic_strain, const double dt,
        const bool update_hist_var) override;

    void setup(const int numgp, const Core::IO::InputParameterContainer& container) override{};

    void pre_evaluate(int gp) override{};

    void update() override{};

    void update_gp_state(int gp) override{};

    void pack_viscoplastic_law(Core::Communication::PackBuffer& data) const override{};

    void unpack_viscoplastic_law(Core::Communication::UnpackBuffer& buffer) override{};
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
