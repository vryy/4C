// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_constraint_springdashpot_manager.hpp"

#include "4C_constraint_springdashpot.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_visualization_parameters.hpp"

FOUR_C_NAMESPACE_OPEN

Constraints::SpringDashpotManager::SpringDashpotManager(
    Global::Problem& problem, std::shared_ptr<Core::FE::Discretization> dis)
    : actdisc_(dis), havespringdashpot_(false)
{
  // get all spring dashpot conditions
  std::vector<const Core::Conditions::Condition*> springdashpots;
  actdisc_->get_condition("RobinSpringDashpot", springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();

  // check if spring dashpots are present
  if (n_conds_)
  {
    // yes! we have spring dashpot conditions!
    havespringdashpot_ = true;

    // new instance of spring dashpot BC with current condition for every spring dashpot condition
    for (int i = 0; i < n_conds_; ++i)
      springs_.push_back(std::make_shared<SpringDashpot>(actdisc_, *springdashpots[i]));

    // set up vtk writer for spring dashpot if requested
    if (problem.io_params().get<bool>("OUTPUT_SPRING"))
    {
      double restart_time = 0.0;
      int restart_step = problem.restart();
      if (restart_step > 0)
      {
        Core::IO::DiscretizationReader reader(
            *actdisc_, problem.input_control_file(), restart_step);
        restart_time = reader.read_double("time");
      }

      const auto use_all_elements = [](const Core::Elements::Element*) { return true; };
      constraint_springdashpot_visualization_writer_ =
          std::make_unique<Core::IO::DiscretizationVisualizationWriterMesh>(actdisc_,
              Core::IO::visualization_parameters_factory(
                  problem.io_params().sublist("RUNTIME VTK OUTPUT"), *problem.output_control_file(),
                  restart_time),
              use_all_elements, actdisc_->name() + "-springdashpot");
    }
  }

  return;
}

void Constraints::SpringDashpotManager::stiffness_and_internal_forces(
    std::shared_ptr<Core::LinAlg::SparseMatrix> stiff,
    std::shared_ptr<Core::LinAlg::Vector<double>> fint,
    std::shared_ptr<Core::LinAlg::Vector<double>> disn, Core::LinAlg::Vector<double>& veln,
    Teuchos::ParameterList parlist)
{
  // evaluate all spring dashpot conditions
  for (int i = 0; i < n_conds_; ++i)
  {
    springs_[i]->reset_newton();
    // get spring type from current condition
    const Constraints::SpringDashpot::RobinSpringDashpotType stype = springs_[i]->get_spring_type();

    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::xyz or
        stype == Constraints::SpringDashpot::RobinSpringDashpotType::refsurfnormal)
      springs_[i]->evaluate_robin(stiff, fint, disn, veln, parlist);
    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::cursurfnormal)
      springs_[i]->evaluate_force_stiff(*stiff, *fint, disn, veln, parlist);
  }

  return;
}

void Constraints::SpringDashpotManager::update()
{
  // update all spring dashpot conditions for each new time step
  for (int i = 0; i < n_conds_; ++i) springs_[i]->update();

  return;
}

void Constraints::SpringDashpotManager::reset_prestress(Core::LinAlg::Vector<double>& dis)
{
  // loop over all spring dashpot conditions and reset them
  for (int i = 0; i < n_conds_; ++i) springs_[i]->reset_prestress(dis);

  return;
}

void Constraints::SpringDashpotManager::output(const double time, const int step)
{
  // only write spring output, if defined in input file
  if (constraint_springdashpot_visualization_writer_ == nullptr) return;

  // row maps for export
  Core::LinAlg::Vector<double> gap(*(actdisc_->node_row_map()), true);
  Core::LinAlg::MultiVector<double> normals(*(actdisc_->node_row_map()), 3, true);
  Core::LinAlg::MultiVector<double> springstress(*(actdisc_->node_row_map()), 3, true);

  // collect outputs from all spring dashpot conditions
  bool found_cursurfnormal = false;
  for (int i = 0; i < n_conds_; ++i)
  {
    springs_[i]->output_gap_normal(gap, normals, springstress);

    // get spring type from current condition
    const Constraints::SpringDashpot::RobinSpringDashpotType stype = springs_[i]->get_spring_type();
    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::cursurfnormal)
      found_cursurfnormal = true;
  }

  // reset time and time step of writer object
  constraint_springdashpot_visualization_writer_->reset();

  // write vectors to output
  if (found_cursurfnormal)
  {
    constraint_springdashpot_visualization_writer_->append_result_data_vector_with_context(
        gap, Core::IO::OutputEntity::node, {"gap"});
    const std::vector<std::optional<std::string>> context(3, "curnormals");
    constraint_springdashpot_visualization_writer_->append_result_data_vector_with_context(
        normals, Core::IO::OutputEntity::node, context);
  }

  // write spring stress
  const std::vector<std::optional<std::string>> context(3, "springstress");
  constraint_springdashpot_visualization_writer_->append_result_data_vector_with_context(
      springstress, Core::IO::OutputEntity::node, context);

  constraint_springdashpot_visualization_writer_->write_to_disk(time, step);
}

void Constraints::SpringDashpotManager::output_restart(
    std::shared_ptr<Core::IO::DiscretizationWriter> output_restart, const double time,
    const int step)
{
  // row maps for export
  std::shared_ptr<Core::LinAlg::Vector<double>> springoffsetprestr =
      std::make_shared<Core::LinAlg::Vector<double>>(*actdisc_->dof_row_map());
  std::shared_ptr<Core::LinAlg::MultiVector<double>> springoffsetprestr_old =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*(actdisc_->node_row_map()), 3, true);

  // collect outputs from all spring dashpot conditions
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const Constraints::SpringDashpot::RobinSpringDashpotType stype = springs_[i]->get_spring_type();

    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::xyz or
        stype == Constraints::SpringDashpot::RobinSpringDashpotType::refsurfnormal)
      springs_[i]->output_prestr_offset(*springoffsetprestr);
    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::cursurfnormal)
      springs_[i]->output_prestr_offset_old(*springoffsetprestr_old);
  }

  // write vector to output for restart
  output_restart->write_vector("springoffsetprestr", springoffsetprestr);
  // write vector to output for restart
  output_restart->write_multi_vector("springoffsetprestr_old", *springoffsetprestr_old);

  // normal output as well
  output(time, step);

  return;
}


/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void Constraints::SpringDashpotManager::read_restart(
    Core::IO::DiscretizationReader& reader, const double& time)
{
  std::shared_ptr<Core::LinAlg::Vector<double>> tempvec =
      std::make_shared<Core::LinAlg::Vector<double>>(*actdisc_->dof_row_map());
  std::shared_ptr<Core::LinAlg::MultiVector<double>> tempvecold =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*(actdisc_->node_row_map()), 3, true);

  reader.read_vector(tempvec, "springoffsetprestr");
  reader.read_multi_vector(tempvecold, "springoffsetprestr_old");

  // loop over all spring dashpot conditions and set restart
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const Constraints::SpringDashpot::RobinSpringDashpotType stype = springs_[i]->get_spring_type();

    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::xyz or
        stype == Constraints::SpringDashpot::RobinSpringDashpotType::refsurfnormal)
      springs_[i]->set_restart(*tempvec);
    if (stype == Constraints::SpringDashpot::RobinSpringDashpotType::cursurfnormal)
      springs_[i]->set_restart_old(*tempvecold);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
