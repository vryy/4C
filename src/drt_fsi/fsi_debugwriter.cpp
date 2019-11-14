/*----------------------------------------------------------------------*/
/*! \file
\brief write debug information for fsi applications
\level 2
\maintainer Matthias Mayr
*----------------------------------------------------------------------*/

#include <sstream>

#include "fsi_debugwriter.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "fsi_monolithic.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::DebugWriter::DebugWriter(Teuchos::RCP<DRT::Discretization> dis) : itnum_(-1)
{
  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  dis_ = discreator->CreateMatchingDiscretizationFromCondition(
      *dis, "FSICoupling", "boundary", "BELE3_3", conditions_to_copy);

  dis_->FillComplete(true, true, true);

  coup_ = Teuchos::rcp(new ADAPTER::Coupling());
  const int ndim = DRT::Problem::Instance()->NDim();
  coup_->SetupCoupling(*dis, *dis_, *DRT::UTILS::ConditionNodeRowMap(*dis, "FSICoupling"),
      *dis_->NodeRowMap(), ndim);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewTimeStep(int step, std::string name)
{
  std::stringstream s;
  s << DRT::Problem::Instance()->OutputControlFile()->FileName();
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new IO::OutputControl(dis_->Comm(),
      "none",                                        // we do not have a problem type
      SHAPEFUNCTION_TYPE::shapefunction_polynomial,  // this is a FE code ... no nurbs
      "debug-output",                                // no input file either
      s.str(),                                       // an output file name is needed
      DRT::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

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
  s << DRT::Problem::Instance()->OutputControlFile()->FileName() << "-" << name_;
  if (name != "") s << "-" << name;
  s << "-step" << step;

  control_ = Teuchos::rcp(new IO::OutputControl(dis_->Comm(),
      "none",                                        // we do not have a problem type
      SHAPEFUNCTION_TYPE::shapefunction_polynomial,  // this is a FE code ... no nurbs
      "debug-output",                                // no input file either
      s.str(),                                       // an output file name is needed
      DRT::Problem::Instance()->NDim(),
      0,     // restart is meaningless here
      1000,  // we never expect to get 1000 iterations
      DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

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
