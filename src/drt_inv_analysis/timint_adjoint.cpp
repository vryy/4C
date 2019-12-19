/*----------------------------------------------------------------------*/
/*! \file
\brief Time integration for a hyperelastic quasi static adjoint problem

\level 3

\maintainer Sebastian Brandstaeter
!*/
#include "timint_adjoint.H"

#include <iostream>
#include "../drt_io/io_pstream.H"
#include "Epetra_SerialDenseVector.h"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_locsys.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_structure.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "Epetra_CrsMatrix.h"


/*----------------------------------------------------------------------*/
/* constructor                                               keh 10/13  */
/*----------------------------------------------------------------------*/
STR::TimIntAdjoint::TimIntAdjoint(Teuchos::RCP<DRT::Discretization> discret)
    : discret_(discret),
      writer_(Teuchos::null),
      solver_(Teuchos::null),
      dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor())),
      dbctoggle_(Teuchos::null),
      dis_(Teuchos::null),
      rhs_(Teuchos::null),
      stiff_(Teuchos::null),
      stiffn_(Teuchos::null),
      zeros_(Teuchos::null),
      stepn_(0),
      dt_(0.0),
      writestep_(0),
      writetime_(0.0),
      isinit_(false)
{
  PrintLogo();

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& io = DRT::Problem::Instance()->IOParams();

  printtoscreen_ = io.get<int>("STDOUTEVRY");

  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  // setup size of the solution vector
  msteps_ = sdyn.get<int>("NUMSTEP");
  dt_ = sdyn.get<double>("TIMESTEP");

  // initialize stiffness matrix
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*(dofrowmap_), 81, true, true));
  stiffn_ = Teuchos::rcp(new LINALG::SparseMatrix(*(dofrowmap_), 81, true, true));

  // initialize zeros from the dofrowmap
  zeros_ = LINALG::CreateVector(*(dofrowmap_), true);

  // initialize solution at stepn
  disdualn_ = LINALG::CreateVector(*(dofrowmap_), true);
  disn_ = LINALG::CreateVector(*(dofrowmap_), true);
  rhsn_ = LINALG::CreateVector(*(dofrowmap_), true);

  dbctoggle_ = LINALG::CreateVector(*(dofrowmap_), true);

  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_, msteps_, true));

  {
    std::vector<DRT::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      locsysman_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_));
    }
  }

  // get the solver associated with structural dynamic section
  CreateLinearSolver(sdyn);

  // Get Dirichlet Map Extractor
  GetDBCMap();

  CreateWriter();
}

void STR::TimIntAdjoint::CreateWriter()
{
  writer_ = Teuchos::rcp(new IO::DiscretizationWriter(discret_));

  // output for the forward problem
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::string prefix = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
  size_t pos = filename.rfind('/');
  size_t pos2 = prefix.rfind('-');
  std::string filenameout = filename.substr(0, pos + 1) + prefix.substr(0, pos2) + "_adjoint" +
                            filename.substr(pos + 1 + prefix.length());

  Teuchos::RCP<IO::OutputControl> controlfile = Teuchos::rcp(new IO::OutputControl(discret_->Comm(),
      DRT::Problem::Instance()->ProblemName(), DRT::Problem::Instance()->SpatialApproximationType(),
      DRT::Problem::Instance()->OutputControlFile()->InputFileName(), filenameout,
      DRT::Problem::Instance()->NDim(), 0,
      DRT::Problem::Instance()->OutputControlFile()->FileSteps(),
      DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

  writer_->SetOutput(controlfile);
  writer_->OverwriteResultFile();
}

/*----------------------------------------------------------------------*/
/* bring primal solution and rhs in                          keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs,
    std::vector<double> mtime, Teuchos::RCP<Epetra_MultiVector> dis, std::vector<double> time)
{
  rhs_ = rhs;
  dis_ = dis;
  time_ = time;
  mtime_ = mtime;

  if ((int)time_.size() != msteps_)
    dserror("setup of the timesteps for the adjoint problem messed up");

  stepn_ = msteps_;
  timen_ = time_[stepn_ - 1];

  if (rhs_ != Teuchos::null && dis_ != Teuchos::null) isinit_ = true;
}

/*----------------------------------------------------------------------*/
/* advance adjoint equation in time                          keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::Integrate()
{
  if (!isinit_) dserror("rhs is not set up properly for the adjoint equation");

  // integration of the adjoints is reverse! For the quasi static case it doesn't matter
  while (stepn_ > 0)
  {
    PrepareStep();
    Solve();
    OutputStep();
    Update();
  }

  isinit_ = false;

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate stiffness at the converged primal state           keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::EvaluateStiff()
{
  // initialize stiffness matrix to zero
  stiff_->Zero();
  stiffn_->Zero();

  // a dummy internal and external force vector
  Teuchos::RCP<Epetra_Vector> fintn = LINALG::CreateVector(*(dofrowmap_), true);
  Teuchos::RCP<Epetra_Vector> fextn = LINALG::CreateVector(*(dofrowmap_), true);

  // set the parameters for the discretization
  Teuchos::ParameterList p;
  p.set("action", "calc_struct_nlnstiff");

  // set time point
  p.set("total time", timen_);
  p.set("delta time", dt_);

  discret_->ClearState();
  discret_->SetState(0, "residual displacement", zeros_);
  discret_->SetState(0, "displacement", disn_);
  discret_->SetState(0, "displacement new", disn_);

  discret_->EvaluateNeumann(p, fextn, stiffn_);
  stiffn_->Complete();

  discret_->Evaluate(p, stiff_, Teuchos::null, fintn, Teuchos::null, Teuchos::null);

  // neumann might bring asymmetrie but we want the adjoint
  // operator to the stiffness matrix, so add transpose
  stiff_->Add(*stiffn_, true, 1.0, 1.0);

  stiff_->Complete();
  discret_->ClearState();
}

Teuchos::RCP<const LINALG::SparseMatrix> STR::TimIntAdjoint::GetLocSysTrafo() const
{
  if (locsysman_ != Teuchos::null) return locsysman_->Trafo();

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* solve the linear system                                   keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::Solve()
{
  // no need to solve for a linear system with rhs of zeros
  if (steprhsn_ == -1)
  {
    disdualn_->Scale(0.0);
    return;
  }

  // transform to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_), rhsn_);

  // Apply Dirichlet
  LINALG::ApplyDirichlettoSystem(
      stiff_, disdualn_, rhsn_, GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

  // Solve
  solver_->Solve(stiff_->EpetraOperator(), disdualn_, rhsn_, true, true, Teuchos::null);

  return;
}

/*----------------------------------------------------------------------*/
/* find step of measurement according to given time          keh 10/14  */
/*----------------------------------------------------------------------*/
int STR::TimIntAdjoint::FindStep(double time)
{
  // find step of the evaluation according to time:
  int step = -1;
  for (int i = 0; i < (int)mtime_.size(); i++)
  {
    double dt = abs(mtime_[i] - time);
    if (dt < 1.0e-10) step = i;
  }

  return step;
}

/*----------------------------------------------------------------------*/
/* prepare step and time                                     keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::PrepareStep()
{
  steprhsn_ = FindStep(timen_);

  // if there is not rhs for this timestep the dual solution is zero anyways
  // so we can skip this
  if (steprhsn_ == -1) return;

  // extract variables needed for this "time step"
  disn_->Update(1.0, *(*dis_)(stepn_ - 1), 0.0);
  rhsn_->Update(-1.0, *(*rhs_)(steprhsn_), 0.0);
  // !!! rhs is the optimality condition differentiated w.r.t the displacement; so it has to be
  // multiplied by -1 to be the correct rhsn_ for the dual problem

  // set internal history needed for this timestep
  SetTimeStepHistory(stepn_);

  // get stiffness matrix according to stepn_
  EvaluateStiff();

  return;
}


/*----------------------------------------------------------------------*/
/* update                                                    keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::Update()
{
  // Update State
  UpdateStepState();

  // Update Step
  stepn_ -= 1;
  timen_ -= dt_;

  if (printtoscreen_ and stepn_ % printtoscreen_ == 0)
  {
    if (discret_->Comm().MyPID() == 0)
      std::cout << "Finalized step " << stepn_ + 1 << " / " << msteps_ << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* update state                                              keh 12/14  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::UpdateStepState()
{
  // Update State
  (*disdual_)(stepn_ - 1)->Update(1.0, *disdualn_, 0.0);

  return;
}

/*----------------------------------------------------------------------*/
/* output                                                    keh 10/14  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::OutputStep()
{
  bool writeown = false;
  if (not writestep_)
  {
    writer_->WriteMesh(0, 0.0);
    writeown = true;
  }

  writestep_ += 1;
  writetime_ += dt_;

  // dont write if the solution is zero
  if (steprhsn_ == -1) return;

  writer_->NewStep(writestep_, writetime_);
  writer_->WriteVector("displacement", disdualn_);
  writer_->WriteElementData(writeown);

  return;
}

/*----------------------------------------------------------------------*/
/* get a linear solver                                       keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::CreateLinearSolver(const Teuchos::ParameterList& sdyn)
{
  // get the solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  solver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
      discret_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  discret_->ComputeNullSpaceIfNecessary(solver_->Params());

  return;
}

/*----------------------------------------------------------------------*/
/* get dirichlet stuff                                       keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::GetDBCMap()
{
  // Get Dirichlet Map
  Teuchos::ParameterList p;
  p.set("total time", timen_);
  discret_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null, dbctoggle_, dbcmaps_);
  zeros_->PutScalar(0.0);  // paranoia

  return;
}

/*----------------------------------------------------------------------*/
/* print logo                                                keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjoint::PrintLogo()
{
  if (discret_->Comm().MyPID() == 0)
  {
    IO::cout << "--------------------------------------------------" << IO::endl;
    IO::cout << "--   Welcome to the adjoint time integration    --" << IO::endl;
    IO::cout << "--------------------------------------------------" << IO::endl;
  }
  return;
}

void STR::TimIntAdjoint::SetTimeStepHistory(int timestep)
{
  Teuchos::ParameterList p;
  p.set("timestep", timestep);
  p.set("action", "calc_struct_recover_istep");
  discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}
