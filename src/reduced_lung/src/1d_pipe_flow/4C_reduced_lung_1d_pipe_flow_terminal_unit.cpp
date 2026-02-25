// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_1d_pipe_flow_terminal_unit.hpp"

#include "4C_reduced_lung_1d_pipe_flow_input.hpp"
#include "4C_utils_exceptions.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN
namespace ReducedLung1DPipe
{
  namespace TerminalUnit
  {

    double evaluate_linear_elastic_pressure(LinearElasticity& linear_elastic_model,
        const TerminalUnitData& data, const double Q, const double dt)
    {
      linear_elastic_model.elastic_pressure_p_el =
          linear_elastic_model.elasticity_E *
          ((data.volume_v + dt * Q) / data.reference_volume_v0 - 1);
      return linear_elastic_model.elastic_pressure_p_el;
    }

    double evaluate_linear_elastic_pressure_gradient(LinearElasticity& linear_elastic_model,
        const TerminalUnitData& data, const double dQdA, const double dt)
    {
      linear_elastic_model.elastic_pressure_grad_dp_el =
          linear_elastic_model.elasticity_E * dt / data.reference_volume_v0 * (dQdA);
      return linear_elastic_model.elastic_pressure_grad_dp_el;
    }

    double evaluate_ogden_hyperelastic_pressure(OgdenHyperelasticity& ogden_hyperelastic_model,
        const double Q, const TerminalUnitData& data, const double dt)
    {
      const double v0_over_v = data.reference_volume_v0 / (data.volume_v + dt * Q);
      ogden_hyperelastic_model.elastic_pressure_p_el =
          ogden_hyperelastic_model.bulk_modulus_kappa /
          ogden_hyperelastic_model.nonlinear_stiffening_beta * v0_over_v *
          (1 - std::pow(v0_over_v, ogden_hyperelastic_model.nonlinear_stiffening_beta));
      return ogden_hyperelastic_model.elastic_pressure_p_el;
    }

    double evaluate_ogden_hyperelastic_pressure_gradient(
        OgdenHyperelasticity& ogden_hyperelastic_model, const TerminalUnitData& data,
        const double Q, const double dQdA, const double dt)
    {
      const double v0_over_v = data.reference_volume_v0 / (data.volume_v + dt * Q);
      ogden_hyperelastic_model.elastic_pressure_grad_dp_el =
          ogden_hyperelastic_model.bulk_modulus_kappa * dt /
          ogden_hyperelastic_model.nonlinear_stiffening_beta *
          (1 - std::pow(v0_over_v, ogden_hyperelastic_model.nonlinear_stiffening_beta) *
                   (1 + ogden_hyperelastic_model.nonlinear_stiffening_beta)) *
          (-v0_over_v / (data.volume_v + dt * Q)) * dQdA;
      return ogden_hyperelastic_model.elastic_pressure_grad_dp_el;
    }

    double evaluate_kelvin_voigt_residual(double area_A, double velocity_u,
        const double reference_area_A0, const double beta, const double Pext,
        const KelvinVoigt& kelvin_voigt_model, const double elastic_pressure_p_el,
        const TerminalUnitData& data)
    {
      // p_1D == p_0D && Q_1D == Q_0D
      const double p_1D = beta * (std::sqrt(area_A) - std::sqrt(reference_area_A0)) + Pext;
      const double flow_Q = area_A * velocity_u;

      const double p_0D = elastic_pressure_p_el +
                          kelvin_voigt_model.viscosity_eta * flow_Q / data.reference_volume_v0;
      return (p_1D - p_0D);
    }

    double evaluate_kelvin_voigt_jacobian(const double area_A, const double dQdA, const double beta,
        const KelvinVoigt& kelvin_voigt_model, const double dp_el_dA, const TerminalUnitData& data)
    {
      // compute (par f) / (par A)
      const double dpdA_1D = 0.5 * beta / std::sqrt(area_A);

      const double dpdA_0D =
          dp_el_dA + kelvin_voigt_model.viscosity_eta / data.reference_volume_v0 * dQdA;
      return dpdA_1D - dpdA_0D;
    }

    double evaluate_four_element_maxwell_residual(const double area_A, const double velocity_u,
        const double reference_area_A0, const double beta, const double Pext,
        const FourElementMaxwell& four_element_maxwell_model, const double elastic_pressure_p_el,
        const double dt, const TerminalUnitData& data)
    {
      // p_1D == p_0D && Q_1D == Q_0D
      const double p_1D = beta * (std::sqrt(area_A) - std::sqrt(reference_area_A0)) + Pext;
      const double flow_Q = area_A * velocity_u;

      const double p_0D = elastic_pressure_p_el +
                          (four_element_maxwell_model.viscosity_eta +
                              (four_element_maxwell_model.elasticity_E_m * dt *
                                  four_element_maxwell_model.viscosity_eta_m) /
                                  (four_element_maxwell_model.elasticity_E_m * dt +
                                      four_element_maxwell_model.viscosity_eta_m)) *
                              flow_Q / data.reference_volume_v0 +
                          four_element_maxwell_model.viscosity_eta_m /
                              (four_element_maxwell_model.elasticity_E_m * dt +
                                  four_element_maxwell_model.viscosity_eta_m) *
                              four_element_maxwell_model.maxwell_pressure_p_m;
      return p_1D - p_0D;
    }

    double evaluate_four_element_maxwell_jacobian(const double area_A, const double dQdA,
        const double beta, const FourElementMaxwell& four_element_maxwell_model,
        const double dp_el_dA, const double dt, const TerminalUnitData& data)
    {
      // compute (par f) / (par A)
      const double dpdA_1D = 0.5 * beta / std::sqrt(area_A);

      const double dpdA_0D = dp_el_dA + (four_element_maxwell_model.viscosity_eta +
                                            (four_element_maxwell_model.elasticity_E_m * dt *
                                                four_element_maxwell_model.viscosity_eta_m) /
                                                (four_element_maxwell_model.elasticity_E_m * dt +
                                                    four_element_maxwell_model.viscosity_eta_m)) *
                                            dQdA / data.reference_volume_v0;
      return dpdA_1D - dpdA_0D;
    }

    ElasticityModel create_elasticity_model(
        const ReducedLung1dPipeFlow::Parameters::TerminalUnits::ElasticityModel& elasticity_model,
        const int global_id)
    {
      if (elasticity_model.elasticity_model_type.at(global_id) ==
          ReducedLung1dPipeFlow::Parameters::TerminalUnits::ElasticityModel::ElasticityModelType::
              Linear)
      {
        return LinearElasticity{.elasticity_E = elasticity_model.linear.elasticity_e.at(global_id),
            .elastic_pressure_p_el = 0.0,
            .elastic_pressure_grad_dp_el = 0.0};
      }
      else if (elasticity_model.elasticity_model_type.at(global_id) ==
               ReducedLung1dPipeFlow::Parameters::TerminalUnits::ElasticityModel::
                   ElasticityModelType::Ogden)
      {
        return OgdenHyperelasticity{
            .bulk_modulus_kappa = elasticity_model.ogden.ogden_parameter_kappa.at(global_id),
            .nonlinear_stiffening_beta = elasticity_model.ogden.ogden_parameter_beta.at(global_id),
            .elastic_pressure_p_el = 0.0,
            .elastic_pressure_grad_dp_el = 0.0};
      }
      FOUR_C_THROW(
          "Unsupported terminal unit elasticity model type for global_id = %d.", global_id);
    }

    RheologicalModel create_rheological_model(
        const ReducedLung1dPipeFlow::Parameters::TerminalUnits::RheologicalModel& rheology_model,
        const int global_id)
    {
      if (rheology_model.rheological_model_type.at(global_id) ==
          ReducedLung1dPipeFlow::Parameters::TerminalUnits::RheologicalModel::RheologicalModelType::
              KelvinVoigt)
      {
        return KelvinVoigt{
            .viscosity_eta = rheology_model.kelvin_voigt.viscosity_kelvin_voigt_eta.at(global_id)};
      }
      else if (rheology_model.rheological_model_type.at(global_id) ==
               ReducedLung1dPipeFlow::Parameters::TerminalUnits::RheologicalModel::
                   RheologicalModelType::FourElementMaxwell)
      {
        return FourElementMaxwell{
            .viscosity_eta =
                rheology_model.four_element_maxwell.viscosity_kelvin_voigt_eta.at(global_id),
            .elasticity_E_m =
                rheology_model.four_element_maxwell.elasticity_maxwell_e_m.at(global_id),
            .viscosity_eta_m =
                rheology_model.four_element_maxwell.viscosity_maxwell_eta_m.at(global_id),
            .maxwell_pressure_p_m = 0.0};
      }
      FOUR_C_THROW(
          "Unsupported terminal unit rheological model type for global_id = %d.", global_id);
    }

    double TerminalUnitModel::evaluate_elastic_pressure(const double area_A, const double beta,
        const double density_rho, const double characteristic_W_outgoing, double dt)
    {
      const double velocity_u =
          characteristic_W_outgoing - 4 * std::pow(area_A, 0.25) * sqrt(0.5 * beta / density_rho);
      double flow_Q = velocity_u * area_A;
      // compute elastic pressure depending on ElasticityModel
      return std::visit(
          [&](auto& elasticity) -> double
          {
            using T = std::decay_t<decltype(elasticity)>;

            if constexpr (std::is_same_v<T, LinearElasticity>)
            {
              return evaluate_linear_elastic_pressure(elasticity, data, flow_Q, dt);
            }
            else if constexpr (std::is_same_v<T, OgdenHyperelasticity>)
            {
              return evaluate_ogden_hyperelastic_pressure(elasticity, flow_Q, data, dt);
            }
          },
          elasticity_model);
    }

    std::pair<double, double> TerminalUnitModel::evaluate_residual_jacobian(double area_A,
        const double reference_area_A0, const double beta, double Pext, const double density_rho,
        const double characteristic_W_outgoing, double dt)
    {
      const double velocity_u =
          characteristic_W_outgoing - 4 * std::pow(area_A, 0.25) * sqrt(0.5 * beta / density_rho);
      double flow_Q = velocity_u * area_A;
      double dQdA = velocity_u - std::pow(area_A, -0.75) * sqrt(0.5 * beta / density_rho) * area_A;

      double p_el = 0.0;
      double dp_eldA = 0.0;
      // compute elastic pressure depending on ElasticityModel
      std::visit(
          [&](auto& elasticity)
          {
            using T = std::decay_t<decltype(elasticity)>;

            if constexpr (std::is_same_v<T, LinearElasticity>)
            {
              p_el = evaluate_linear_elastic_pressure(elasticity, data, flow_Q, dt);
              dp_eldA = evaluate_linear_elastic_pressure_gradient(elasticity, data, dQdA, dt);
            }
            else if constexpr (std::is_same_v<T, OgdenHyperelasticity>)
            {
              p_el = evaluate_ogden_hyperelastic_pressure(elasticity, flow_Q, data, dt);
              dp_eldA =
                  evaluate_ogden_hyperelastic_pressure_gradient(elasticity, data, flow_Q, dQdA, dt);
            }
          },
          elasticity_model);

      // compute residual and jacobian depending on rheological model
      return std::visit(
          [&](auto& rheology) -> std::pair<double, double>
          {
            using T = std::decay_t<decltype(rheology)>;

            if constexpr (std::is_same_v<T, KelvinVoigt>)
            {
              double residual = evaluate_kelvin_voigt_residual(
                  area_A, velocity_u, reference_area_A0, beta, Pext, rheology, p_el, data);
              double jacobian =
                  evaluate_kelvin_voigt_jacobian(area_A, dQdA, beta, rheology, dp_eldA, data);
              return std::make_pair(residual, jacobian);
            }
            else if constexpr (std::is_same_v<T, FourElementMaxwell>)
            {
              double residual = evaluate_four_element_maxwell_residual(
                  area_A, velocity_u, reference_area_A0, beta, Pext, rheology, p_el, dt, data);
              double jacobian = evaluate_four_element_maxwell_jacobian(
                  area_A, dQdA, beta, rheology, dp_eldA, dt, data);
              return std::make_pair(residual, jacobian);
            }
            else
            {
              FOUR_C_ASSERT(false, "Wrongly defined rheological model.");
            }
          },
          rheological_model);
    }

    void TerminalUnitModel::update_terminal_unit_data(const double flow_Q, const double dt,
        const double area_A, const double beta, const double density_rho,
        const double characteristic_W_outgoing)
    {
      data.volume_v += dt * flow_Q;

      // update p_m in FourElementMaxwellElement
      auto* four_element_maxwell = std::get_if<FourElementMaxwell>(&rheological_model);
      if (four_element_maxwell)
      {
        double p_m = four_element_maxwell->maxwell_pressure_p_m;
        double p_el =
            evaluate_elastic_pressure(area_A, beta, density_rho, characteristic_W_outgoing, dt);
        four_element_maxwell->maxwell_pressure_p_m =
            four_element_maxwell->viscosity_eta_m /
                (four_element_maxwell->elasticity_E_m * dt +
                    four_element_maxwell->viscosity_eta_m) *
                p_m +
            four_element_maxwell->elasticity_E_m * dt /
                (four_element_maxwell->elasticity_E_m * dt +
                    four_element_maxwell->viscosity_eta_m) *
                p_el;
      }
    }

  }  // namespace TerminalUnit
}  // namespace ReducedLung1DPipe
FOUR_C_NAMESPACE_CLOSE