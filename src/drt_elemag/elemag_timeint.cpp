/*----------------------------------------------------------------------*/
/*!
\file elemag_timeint.cpp

\brief Base class functions for time integration of electromagnetics

<pre>
\level 3

\maintainer Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "elemag_timeint.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"


#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
ELEMAG::ElemagTimeInt::ElemagTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
    : discret_(actdis),
      solver_(solver),
      params_(params),
      output_(output),
      elemagdyna_(DRT::INPUT::IntegralValue<INPAR::ELEMAG::DynamicType>(*params_, "TIMEINT")),
      myrank_(actdis->Comm().MyPID()),
      time_(0.0),
      step_(0),
      restart_(params_->get<int>("restart")),
      maxtime_(params_->get<double>("MAXTIME")),
      stepmax_(params_->get<int>("NUMSTEP")),
      uprestart_(params_->get<int>("RESTARTEVRY", -1)),
      upres_(params_->get<int>("RESULTSEVRY", -1)),
      numdim_(DRT::Problem::Instance()->NDim()),
      dtp_(params_->get<double>("TIMESTEP")),
      dtele_(0.0),
      dtsolve_(0.0),
      calcerr_(false)
{
  // constructor supposed to be empty!

}  // ElemagTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
ELEMAG::ElemagTimeInt::~ElemagTimeInt() {}


/*----------------------------------------------------------------------*
 |  initialization routine (public)                    gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Init()
{
  // get dof row map
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // check time-step length
  if (dtp_ <= 0.0) dserror("Zero or negative time-step length!");

  // create vector of zeros to be used for enforcing zero Dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap, true);

  /*dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);
  }*/

  // create system matrix and set to zero
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 108, false, true));
  sysmat_->Zero();

  // create residual vector
  residual_ = LINALG::CreateVector(*dofrowmap, true);

  // write mesh
  output_->WriteMesh(0, 0.0);

  return;
}  // Init


/*----------------------------------------------------------------------*
 |  Print information to screen (public)               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::PrintInformationToScreen() { return; }  //  PrintInformationToScreen


/*----------------------------------------------------------------------*
 |  Time integration (public)                          gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Integrate()
{
  // time measurement: integration
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::Integrate");

  // output of initial field
  Output();

  // call elements to calculate system matrix/rhs and assemble
  AssembleMatAndRHS();

  // apply Dirichlet boundary conditions to system of equations
  ApplyDirichletToSystem();

  // time loop
  while (step_ < stepmax_ and time_ < maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    // solve
    Solve();

    // output of solution
    Output();

  }  // while (step_<stepmax_ and time_<maxtime_)

  return;
}  // Integrate


/*----------------------------------------------------------------------*
 |  Set initial field by given function (public)       gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::SetInitialField(int startfuncno) { return; }  // SetInitialField


/*----------------------------------------------------------------------*
 |  Assemble matrix and right-hand side (public)       gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::AssembleMatAndRHS() { return; }  // AssembleMatAndRHS


/*----------------------------------------------------------------------*
 |  Apply Dirichlet b.c. to system (public)            gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ApplyDirichletToSystem() { return; }  // ApplyDirichletToSystem


/*----------------------------------------------------------------------*
 |  Solve system for trace and then interior field     gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Solve() { return; }  // Solve


/*----------------------------------------------------------------------*
 |  Output to screen (public)                          gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::OutputToScreen() { return; }  // OutputToScreen


/*----------------------------------------------------------------------*
 |  Output (public)                                    gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Output() { return; }  // Output


/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                     gravemeier 06/17 |
 *----------------------------------------------------------------------*/
/*void ELEMAG::ElemagTimeInt::WriteRestart()
{
  if (myrank_ == 0) std::cout<<"======= Restart written in step " << step_ << std::endl;

  return;
}*/ // WriteRestart


/*----------------------------------------------------------------------*
 |  ReadRestart (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
/*void ELEMAG::ElemagTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  return;
}*/ // ReadRestart



/*----------------------------------------------------------------------*
 |  Create test field (public)                         gravemeier 06/17 |
 *----------------------------------------------------------------------*/
/*Teuchos::RCP<DRT::ResultTest> ELEMAG::ElemagTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new ElemagResultTest(*this));
}*/ // CreateFieldTest
