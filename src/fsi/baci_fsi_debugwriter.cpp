/*----------------------------------------------------------------------*/
/*! \file
\brief write debug information for fsi applications
\level 2
*----------------------------------------------------------------------*/

#include "baci_fsi_debugwriter.H"

#include "baci_adapter_ale_fsi.H"
#include "baci_adapter_fld_fluid_fsi.H"
#include "baci_adapter_str_fsiwrapper.H"
#include "baci_coupling_adapter.H"
#include "baci_fsi_monolithic.H"
#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_lib_discret.H"
#include "baci_lib_utils_createdis.H"

#include <sstream>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::DebugWriter::DebugWriter(Teuchos::RCP<DRT::Discretization> dis) : itnum_(-1)
{
  std::vector<std::string> conditions_to_copy = {"FSICoupling"};
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  dis_ = discreator->CreateMatchingDiscretizationFromCondition(
      *dis, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  dis_->FillComplete(true, true, true);

  coup_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  const int ndim = GLOBAL::Problem::Instance()->NDim();
  coup_->SetupCoupling(*dis, *dis_, *DRT::UTILS::ConditionNodeRowMap(*dis, "FSICoupling"),
      *dis_->NodeRowMap(), ndim);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewTimeStep(int step, std::string name)
{
  std::stringstream s;
  s << GLOBAL::Problem::Instance()->OutputControlFile()->FileName();
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new IO::OutputControl(dis_->Comm(),
      "none",                                   // we do not have a problem type
      CORE::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      GLOBAL::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      INPUT::IntegralValue<int>(GLOBAL::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

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
    Teuchos::RCP<DRT::Discretization> dis, const std::string& name)
    : dis_(dis), name_(name), itnum_(-1)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SimpleDebugWriter::NewLinearSystem(int step, std::string name)
{
  std::stringstream s;
  s << GLOBAL::Problem::Instance()->OutputControlFile()->FileName() << "-" << name_;
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new IO::OutputControl(dis_->Comm(),
      "none",                                   // we do not have a problem type
      CORE::FE::ShapeFunctionType::polynomial,  // this is a FE code ... no nurbs
      "debug-output",                           // no input file either
      s.str(),                                  // an output file name is needed
      GLOBAL::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      INPUT::IntegralValue<int>(GLOBAL::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

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
      new SimpleDebugWriter(algorithm_.StructureField()->Discretization(), "structure"));
  fluid_writer_ =
      Teuchos::rcp(new SimpleDebugWriter(algorithm_.FluidField()->Discretization(), "fluid"));
  ale_writer_ = Teuchos::rcp(
      new SimpleDebugWriter(algorithm_.AleField()->WriteAccessDiscretization(), "ale"));
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

  algorithm_.ExtractFieldVectors(v, sx, fx, ax);

  Epetra_Vector s(*sx);
  Epetra_Vector f(*fx);
  Epetra_Vector a(*ax);

  struct_writer_->WriteVector(name, s);
  fluid_writer_->WriteVector(name, f);
  ale_writer_->WriteVector(name, a);
}

BACI_NAMESPACE_CLOSE
