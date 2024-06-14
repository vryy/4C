/*----------------------------------------------------------------------*/
/*! \file
\brief write debug information for fsi applications
\level 2
*----------------------------------------------------------------------*/

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
FSI::UTILS::DebugWriter::DebugWriter(Teuchos::RCP<Core::FE::Discretization> dis) : itnum_(-1)
{
  std::vector<std::string> conditions_to_copy = {"FSICoupling"};
  Teuchos::RCP<Core::FE::DiscretizationCreatorBase> discreator =
      Teuchos::rcp(new Core::FE::DiscretizationCreatorBase());
  dis_ = discreator->create_matching_discretization_from_condition(
      *dis, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  dis_->fill_complete(true, true, true);

  coup_ = Teuchos::rcp(new Core::Adapter::Coupling());
  const int ndim = Global::Problem::Instance()->NDim();
  coup_->setup_coupling(*dis, *dis_, *Core::Conditions::ConditionNodeRowMap(*dis, "FSICoupling"),
      *dis_->NodeRowMap(), ndim);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewTimeStep(int step, std::string name)
{
  std::stringstream s;
  s << Global::Problem::Instance()->OutputControlFile()->file_name();
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new Core::IO::OutputControl(dis_->Comm(),
      "none",                                   // we do not have a problem type
      Core::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      Global::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      Core::UTILS::IntegralValue<bool>(Global::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

  writer_ = dis_->Writer();
  writer_->SetOutput(control_);
  itnum_ = 0;
  writer_->WriteMesh(0, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewIteration()
{
  writer_->NewStep(itnum_, itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::WriteVector(const std::string& name, const Epetra_Vector& v)
{
  writer_->WriteVector(name, coup_->MasterToSlave(Teuchos::rcp(&v, false)));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::SimpleDebugWriter::SimpleDebugWriter(
    Teuchos::RCP<Core::FE::Discretization> dis, const std::string& name)
    : dis_(dis), name_(name), itnum_(-1)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SimpleDebugWriter::NewLinearSystem(int step, std::string name)
{
  std::stringstream s;
  s << Global::Problem::Instance()->OutputControlFile()->file_name() << "-" << name_;
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new Core::IO::OutputControl(dis_->Comm(),
      "none",                                   // we do not have a problem type
      Core::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      Global::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      Core::UTILS::IntegralValue<bool>(Global::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

  writer_ = dis_->Writer();
  writer_->SetOutput(control_);
  itnum_ = 0;
  writer_->WriteMesh(0, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SimpleDebugWriter::NewIteration()
{
  writer_->NewStep(itnum_, itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SimpleDebugWriter::WriteVector(const std::string& name, Epetra_Vector& v)
{
  writer_->WriteVector(name, Teuchos::rcp(&v, false));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::MonolithicDebugWriter::MonolithicDebugWriter(Monolithic& algorithm)
    : algorithm_(algorithm), counter_(0)
{
  struct_writer_ = Teuchos::rcp(
      new SimpleDebugWriter(algorithm_.structure_field()->discretization(), "structure"));
  fluid_writer_ =
      Teuchos::rcp(new SimpleDebugWriter(algorithm_.fluid_field()->discretization(), "fluid"));
  ale_writer_ = Teuchos::rcp(
      new SimpleDebugWriter(algorithm_.ale_field()->write_access_discretization(), "ale"));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::MonolithicDebugWriter::NewLinearSystem()
{
  counter_ += 1;
  struct_writer_->NewLinearSystem(counter_);
  fluid_writer_->NewLinearSystem(counter_);
  ale_writer_->NewLinearSystem(counter_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::MonolithicDebugWriter::NewIteration()
{
  struct_writer_->NewIteration();
  fluid_writer_->NewIteration();
  ale_writer_->NewIteration();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::MonolithicDebugWriter::WriteVector(
    const std::string& name, const Teuchos::RCP<Epetra_Vector>& v)
{
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  algorithm_.extract_field_vectors(v, sx, fx, ax);

  Epetra_Vector s(*sx);
  Epetra_Vector f(*fx);
  Epetra_Vector a(*ax);

  struct_writer_->WriteVector(name, s);
  fluid_writer_->WriteVector(name, f);
  ale_writer_->WriteVector(name, a);
}

FOUR_C_NAMESPACE_CLOSE
