#ifdef CCADISCRET

#include <sstream>

#include "fsi_debugwriter.H"

#include "../drt_adapter/adapter_utils.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_globalproblem.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::DebugWriter::DebugWriter(Teuchos::RCP<DRT::Discretization> dis)
{
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  dis_ = DRT::UTILS::CreateDiscretizationFromCondition(dis,"FSICoupling","boundary","BELE3",conditions_to_copy);
  dis_->FillComplete();

  dis_->Print(cout);

  coup_.SetupCoupling(*dis,
                      *dis_,
                      *ADAPTER::UTILS::ConditionNodeMap(*dis,"FSICoupling"),
                      *dis_->NodeRowMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewTimeStep(int step, std::string name)
{
  std::stringstream s;
  s << DRT::Problem::Instance()->OutputControlFile()->FileName();
  if (name!="")
    s << "-" << name;
  s << "-step"
    << step;

  control_ = Teuchos::rcp(
    new IO::OutputControl(
      dis_->Comm(),
      "none",                   // we do not have a problem type
      "Polynomial",             // this is a FE code ... no nurbs
      "debug-output",           // no input file either
      s.str(),                  // an output file name is needed
      genprob.ndim,
      0,                        // restart is meaningless here
      1000));                   // we never expect to get 1000 iterations

  writer_ = Teuchos::rcp(new IO::DiscretizationWriter(dis_,control_));
  itnum_ = 0;
  writer_->WriteMesh(0,0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::NewIteration()
{
  writer_->NewStep(itnum_,itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DebugWriter::WriteVector(const std::string& name, const Epetra_Vector& v)
{
  writer_->WriteVector(name,coup_.MasterToSlave(Teuchos::rcp(&v,false)));
}


#endif
