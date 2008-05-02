#ifdef CCADISCRET

#include <sstream>

#include "fsi_debugwriter.H"

#include "../drt_adapter/adapter_create_boundary.H"
#include "../drt_adapter/adapter_utils.H"

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
FSI::DebugWriter::DebugWriter(Teuchos::RCP<DRT::Discretization> dis)
{
  dis_ = CreateDiscretizationFromCondition(dis,"FSICoupling","boundary","BELE3");
  dis_->FillComplete();

  coup_.SetupCoupling(*dis,
                      *dis_,
                      *ADAPTER::UTILS::ConditionNodeMap(*dis,"FSICoupling"),
                      *dis_->NodeRowMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DebugWriter::NewTimeStep(int step)
{
  std::stringstream s;
  s << DRT::Problem::Instance()->OutputControlFile()->FileName()
    << "-step"
    << step;

  control_ = Teuchos::rcp(
    new IO::OutputControl(
      dis_->Comm(),
      "none",                   // we do not have a problem type
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
void FSI::DebugWriter::NewIteration()
{
  writer_->NewStep(itnum_,itnum_);
  itnum_ += 1;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DebugWriter::WriteVector(const std::string& name, const Epetra_Vector& v)
{
  writer_->WriteVector(name,coup_.MasterToSlave(Teuchos::rcp(&v,false)));
}


#endif
