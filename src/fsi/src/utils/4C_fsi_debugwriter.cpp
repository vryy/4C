// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_debugwriter.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Utils::DebugWriter::DebugWriter(std::shared_ptr<Core::FE::Discretization> dis) : itnum_(-1)
{
  std::vector<std::string> conditions_to_copy = {"FSICoupling"};
  std::shared_ptr<Core::FE::DiscretizationCreatorBase> discreator =
      std::make_shared<Core::FE::DiscretizationCreatorBase>();
  dis_ = discreator->create_matching_discretization_from_condition(
      *dis, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  dis_->fill_complete(true, true, true);

  coup_ = std::make_shared<Coupling::Adapter::Coupling>();
  const int ndim = Global::Problem::instance()->n_dim();
  coup_->setup_coupling(*dis, *dis_, *Core::Conditions::condition_node_row_map(*dis, "FSICoupling"),
      *dis_->node_row_map(), ndim);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::DebugWriter::new_time_step(int step, std::string name)
{
  std::stringstream s;
  s << Global::Problem::instance()->output_control_file()->file_name();
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = std::make_shared<Core::IO::OutputControl>(dis_->get_comm(),
      "none",                                   // we do not have a problem type
      Core::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      Global::Problem::instance()->n_dim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      Global::Problem::instance()->io_params().get<bool>("OUTPUT_BIN"));

  writer_ = dis_->writer();
  writer_->set_output(control_);
  itnum_ = 0;
  writer_->write_mesh(0, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::DebugWriter::new_iteration()
{
  writer_->new_step(itnum_, itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::DebugWriter::write_vector(
    const std::string& name, const Core::LinAlg::Vector<double>& v)
{
  writer_->write_vector(name, coup_->master_to_slave(v));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Utils::SimpleDebugWriter::SimpleDebugWriter(
    std::shared_ptr<Core::FE::Discretization> dis, const std::string& name)
    : dis_(dis), name_(name), itnum_(-1)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::SimpleDebugWriter::new_linear_system(int step, std::string name)
{
  std::stringstream s;
  s << Global::Problem::instance()->output_control_file()->file_name() << "-" << name_;
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = std::make_shared<Core::IO::OutputControl>(dis_->get_comm(),
      "none",                                   // we do not have a problem type
      Core::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      Global::Problem::instance()->n_dim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      Global::Problem::instance()->io_params().get<bool>("OUTPUT_BIN"));

  writer_ = dis_->writer();
  writer_->set_output(control_);
  itnum_ = 0;
  writer_->write_mesh(0, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::SimpleDebugWriter::new_iteration()
{
  writer_->new_step(itnum_, itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::SimpleDebugWriter::write_vector(
    const std::string& name, Core::LinAlg::Vector<double>& v)
{
  writer_->write_vector(name, Core::Utils::shared_ptr_from_ref(v));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Utils::MonolithicDebugWriter::MonolithicDebugWriter(Monolithic& algorithm)
    : algorithm_(algorithm), counter_(0)
{
  struct_writer_ = std::make_shared<SimpleDebugWriter>(
      algorithm_.structure_field()->discretization(), "structure");
  fluid_writer_ =
      std::make_shared<SimpleDebugWriter>(algorithm_.fluid_field()->discretization(), "fluid");
  ale_writer_ = std::make_shared<SimpleDebugWriter>(
      algorithm_.ale_field()->write_access_discretization(), "ale");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::MonolithicDebugWriter::new_linear_system()
{
  counter_ += 1;
  struct_writer_->new_linear_system(counter_);
  fluid_writer_->new_linear_system(counter_);
  ale_writer_->new_linear_system(counter_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::MonolithicDebugWriter::new_iteration()
{
  struct_writer_->new_iteration();
  fluid_writer_->new_iteration();
  ale_writer_->new_iteration();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Utils::MonolithicDebugWriter::write_vector(
    const std::string& name, const std::shared_ptr<Core::LinAlg::Vector<double>>& v)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> sx;
  std::shared_ptr<const Core::LinAlg::Vector<double>> fx;
  std::shared_ptr<const Core::LinAlg::Vector<double>> ax;

  algorithm_.extract_field_vectors(v, sx, fx, ax);

  Core::LinAlg::Vector<double> s(*sx);
  Core::LinAlg::Vector<double> f(*fx);
  Core::LinAlg::Vector<double> a(*ax);

  struct_writer_->write_vector(name, s);
  fluid_writer_->write_vector(name, f);
  ale_writer_->write_vector(name, a);
}

FOUR_C_NAMESPACE_CLOSE
