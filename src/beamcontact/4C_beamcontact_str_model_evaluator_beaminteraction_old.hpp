// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMCONTACT_STR_MODEL_EVALUATOR_BEAMINTERACTION_OLD_HPP
#define FOUR_C_BEAMCONTACT_STR_MODEL_EVALUATOR_BEAMINTERACTION_OLD_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace CONTACT
{
  class Beam3cmanager;
}  // namespace CONTACT

namespace Solid
{
  namespace ModelEvaluator
  {
    class BeamInteractionOld : public Generic
    {
     public:
      //! constructor
      BeamInteractionOld();

      void setup() override;

      //! derived
      Inpar::Solid::ModelType type() const override
      {
        return Inpar::Solid::model_beam_interaction_old;
      }

      //! derived
      void reset(const Core::LinAlg::Vector<double>& x) override;

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override { return; };

      //! derived
      void post_evaluate() override{/* currently unused */};

      //! derived
      bool assemble_force(Core::LinAlg::Vector<double>& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void predict(const Inpar::Solid::PredEnum& pred_type) override { return; };

      //! derived
      void run_pre_compute_x(const Core::LinAlg::Vector<double>& xold,
          Core::LinAlg::Vector<double>& dir_mutable, const NOX::Nln::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void run_post_compute_x(const Core::LinAlg::Vector<double>& xold,
          const Core::LinAlg::Vector<double>& dir,
          const Core::LinAlg::Vector<double>& xnew) override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void update_step_state(const double& timefac_n) override;

      //! derived
      void update_step_element() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void determine_energy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void reset_step_state() override;

      //! [derived]
      void post_output() override;

      //! derived
      std::shared_ptr<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      std::shared_ptr<const Core::LinAlg::Vector<double>> get_current_solution_ptr() const override;

      //! derived
      std::shared_ptr<const Core::LinAlg::Vector<double>> get_last_time_step_solution_ptr()
          const override;

     private:
      //! structural displacement at \f$t_{n+1}\f$
      std::shared_ptr<const Core::LinAlg::Vector<double>> disnp_ptr_;

      //! stiffness contributions from beam interactions
      std::shared_ptr<Core::LinAlg::SparseMatrix> stiff_beaminteract_ptr_;

      //! force contributions from beam interaction at \f$t_{n+1}\f$
      std::shared_ptr<Core::LinAlg::Vector<double>> f_beaminteract_np_ptr_;

      //! beam contact manager
      std::shared_ptr<CONTACT::Beam3cmanager> beamcman_;
    };

  }  // namespace ModelEvaluator
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
