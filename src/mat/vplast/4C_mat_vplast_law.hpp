// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VPLAST_LAW_HPP
#define FOUR_C_MAT_VPLAST_LAW_HPP

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
  namespace Viscoplastic
  {

    namespace PAR
    {

    };

    /*----------------------------------------------------------------------*/
    /*! \class Law
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
    class Law
    {
     public:
      /// construct viscoplastic laws with specific material params
      explicit Law(Core::Mat::PAR::Parameter* params);
      /// construct empty viscoplastic law
      Law();

      /// virtual destructor
      virtual ~Law() = default;

      /*!
       * @brief create object by input parameter ID
       *
       * @param[in] matnum  material ID
       * @return pointer to material that is defined by material ID
       */
      static std::shared_ptr<Law> factory(int matnum);

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
       * @brief Evaluate the equivalent plastic strain rate for a given equivalent stress and a
       * given plastic strain
       *
       * @param[in] equiv_stress Equivalent stress
       * @param[in] equiv_plastic_strain Equivalent plastic strain
       * @param[in] dt Time step size (used solely for overflow checking purposes)
       * @return Equivalent plastic strain rate
       */
      virtual double evaluate_plastic_strain_rate(const double equiv_stress,
          const double equiv_plastic_strain, const double dt,
          const bool update_hist_var = true) = 0;

      /*!
       * @brief Evaluate the derivatives of the equivalent plastic strain rate for a given
       * equivalent stress and a given plastic strain
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
       * @brief Update the history variables for a specific GP after a converged substep
       * (substepping procedure)
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
  }  // namespace Viscoplastic


}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
