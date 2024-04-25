/*----------------------------------------------------------------------*/
/*! \file

\brief Class controlling constraints and containing the necessary data, code originally by Thomas
Kloeppel


\level 2


*----------------------------------------------------------------------*/


#include "4C_constraint_manager.hpp"

#include "4C_constraint.hpp"
#include "4C_constraint_dofset.hpp"
#include "4C_constraint_monitor.hpp"
#include "4C_constraint_multipointconstraint2.hpp"
#include "4C_constraint_multipointconstraint3.hpp"
#include "4C_constraint_multipointconstraint3penalty.hpp"
#include "4C_constraint_penalty.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONSTRAINTS::ConstrManager::ConstrManager()
    : offset_id_(-1),
      max_constr_id_(0),
      num_constr_id_(-1),
      num_monitor_id_(-1),
      min_monitor_id_(-1),
      haveconstraint_(false),
      havelagrconstr_(false),
      havepenaconstr_(false),
      havemonitor_(false),
      uzawaparam_(0.0),
      issetup_(false),
      isinit_(false)

{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::Init(
    Teuchos::RCP<DRT::Discretization> discr, const Teuchos::ParameterList& params)
{
  SetIsSetup(false);

  // set pointer to discretization
  actdisc_ = discr;

  //----------------------------------------------------------------------------
  //---------------------------------------------------------Constraint Conditions!

  // constructors of constraints increment number of constraints defined and the minimum
  // ConditionID read so far.
  num_constr_id_ = 0;
  offset_id_ = 10000;
  // Check, what kind of constraining boundary conditions there are
  volconstr3d_ =
      Teuchos::rcp(new Constraint(actdisc_, "VolumeConstraint_3D", offset_id_, max_constr_id_));
  areaconstr3d_ =
      Teuchos::rcp(new Constraint(actdisc_, "AreaConstraint_3D", offset_id_, max_constr_id_));
  areaconstr2d_ =
      Teuchos::rcp(new Constraint(actdisc_, "AreaConstraint_2D", offset_id_, max_constr_id_));
  mpconline2d_ =
      Teuchos::rcp(new MPConstraint2(actdisc_, "MPC_NodeOnLine_2D", offset_id_, max_constr_id_));
  mpconplane3d_ =
      Teuchos::rcp(new MPConstraint3(actdisc_, "MPC_NodeOnPlane_3D", offset_id_, max_constr_id_));
  mpcnormcomp3d_ = Teuchos::rcp(
      new MPConstraint3(actdisc_, "MPC_NormalComponent_3D", offset_id_, max_constr_id_));

  volconstr3dpen_ = Teuchos::rcp(new ConstraintPenalty(actdisc_, "VolumeConstraint_3D_Pen"));
  areaconstr3dpen_ = Teuchos::rcp(new ConstraintPenalty(actdisc_, "AreaConstraint_3D_Pen"));
  mpcnormcomp3dpen_ =
      Teuchos::rcp(new MPConstraint3Penalty(actdisc_, "MPC_NormalComponent_3D_Pen"));

  havepenaconstr_ = (mpcnormcomp3dpen_->HaveConstraint()) or (volconstr3dpen_->HaveConstraint()) or
                    (areaconstr3dpen_->HaveConstraint());

  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  havelagrconstr_ = (areaconstr3d_->HaveConstraint()) or (volconstr3d_->HaveConstraint()) or
                    (areaconstr2d_->HaveConstraint()) or (mpconplane3d_->HaveConstraint()) or
                    (mpcnormcomp3d_->HaveConstraint()) or (mpconline2d_->HaveConstraint());
  haveconstraint_ = havepenaconstr_ or havelagrconstr_;


  SetIsInit(true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::Setup(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::ParameterList params)
{
  CheckIsInit();

  if (haveconstraint_)
  {
    num_constr_id_ = std::max(max_constr_id_ - offset_id_ + 1, 0);
    constrdofset_ = Teuchos::rcp(new ConstraintDofSet());
    constrdofset_->AssignDegreesOfFreedom(actdisc_, num_constr_id_, 0);
    offset_id_ -= constrdofset_->FirstGID();
    Teuchos::ParameterList p;
    uzawaparam_ = params.get<double>("uzawa parameter", 1);
    double time = params.get<double>("total time", 0.0);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    // initialize constrMatrix
    constr_matrix_ =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(*dofrowmap, num_constr_id_, false, true));
    // build Epetra_Map used as domainmap for constrMatrix and rowmap for result vectors
    constrmap_ = Teuchos::rcp(new Epetra_Map(*(constrdofset_->DofRowMap())));
    // build an all reduced version of the constraintmap, since sometimes all processors
    // have to know all values of the constraints and Lagrange multipliers
    redconstrmap_ = CORE::LINALG::AllreduceEMap(*constrmap_);
    // importer
    conimpo_ = Teuchos::rcp(new Epetra_Export(*redconstrmap_, *constrmap_));
    // sum up initial values
    refbasevalues_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
    Teuchos::RCP<Epetra_Vector> refbaseredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
    // Compute initial values and assemble them to the completely redundant vector
    // We will always use the third systemvector for this purpose
    p.set("OffsetID", offset_id_);
    p.set("total time", time);
    actdisc_->SetState("displacement", disp);
    volconstr3d_->Initialize(p, refbaseredundant);
    areaconstr3d_->Initialize(p, refbaseredundant);
    areaconstr2d_->Initialize(p, refbaseredundant);
    volconstr3dpen_->Initialize(p);
    areaconstr3dpen_->Initialize(p);

    mpconline2d_->SetConstrState("displacement", disp);
    mpconline2d_->Initialize(p, refbaseredundant);
    mpconplane3d_->SetConstrState("displacement", disp);
    mpconplane3d_->Initialize(p, refbaseredundant);
    mpcnormcomp3d_->SetConstrState("displacement", disp);
    mpcnormcomp3d_->Initialize(p, refbaseredundant);
    mpcnormcomp3dpen_->SetConstrState("displacement", disp);
    mpcnormcomp3dpen_->Initialize(p);

    // Export redundant vector into distributed one
    refbasevalues_->Export(*refbaseredundant, *conimpo_, Add);

    // Initialize Lagrange Multipliers, reference values and errors
    actdisc_->ClearState();
    referencevalues_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
    actvalues_ = Teuchos::rcp(new Epetra_Vector(*constrmap_, true));
    constrainterr_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
    lagr_mult_vec_ = Teuchos::rcp(new Epetra_Vector(*constrmap_, true));
    lagr_mult_vec_old_ = Teuchos::rcp(new Epetra_Vector(*constrmap_, true));
    fact_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
  }
  //----------------------------------------------------------------------------
  //---------------------------------------------------------Monitor Conditions!
  actdisc_->SetState("displacement", disp);
  min_monitor_id_ = 10000;
  int maxMonitorID = 0;
  volmonitor3d_ =
      Teuchos::rcp(new Monitor(actdisc_, "VolumeMonitor_3D", min_monitor_id_, maxMonitorID));
  areamonitor3d_ =
      Teuchos::rcp(new Monitor(actdisc_, "AreaMonitor_3D", min_monitor_id_, maxMonitorID));
  areamonitor2d_ =
      Teuchos::rcp(new Monitor(actdisc_, "AreaMonitor_2D", min_monitor_id_, maxMonitorID));
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  num_monitor_id_ = std::max(maxMonitorID - min_monitor_id_ + 1, 0);
  havemonitor_ = (areamonitor3d_->HaveMonitor()) || (volmonitor3d_->HaveMonitor()) ||
                 (areamonitor2d_->HaveMonitor());
  if (havemonitor_)
  {
    Teuchos::ParameterList p1;
    // monitor values are only stored on processor zero since they are only needed for output
    int nummyele = 0;
    if (!actdisc_->Comm().MyPID())
    {
      nummyele = num_monitor_id_;
    }
    // initialize maps and importer
    monitormap_ = Teuchos::rcp(new Epetra_Map(num_monitor_id_, nummyele, 0, actdisc_->Comm()));
    redmonmap_ = CORE::LINALG::AllreduceEMap(*monitormap_);
    monimpo_ = Teuchos::rcp(new Epetra_Export(*redmonmap_, *monitormap_));
    monitorvalues_ = Teuchos::rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_ = Teuchos::rcp(new Epetra_Vector(*monitormap_));

    Teuchos::RCP<Epetra_Vector> initialmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
    p1.set("OffsetID", min_monitor_id_);
    volmonitor3d_->Evaluate(p1, initialmonredundant);
    areamonitor3d_->Evaluate(p1, initialmonredundant);
    areamonitor2d_->Evaluate(p1, initialmonredundant);

    // Export redundant vector into distributed one
    initialmonvalues_->Export(*initialmonredundant, *monimpo_, Add);
    monitortypes_ = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
    BuildMoniType();
  }

  SetIsSetup(true);
}


/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::EvaluateForceStiff(const double time,
    Teuchos::RCP<const Epetra_Vector> displast, Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> fint, Teuchos::RCP<CORE::LINALG::SparseOperator> stiff,
    Teuchos::ParameterList scalelist)
{
  CheckIsInit();
  CheckIsSetup();

  double scStiff = scalelist.get("scaleStiffEntries", 1.0);
  double scConMat = scalelist.get("scaleConstrMat", 1.0);

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  std::vector<DRT::Condition*> constrcond(0);
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
  constr_matrix_->Reset();  //=Teuchos::rcp(new
                            // CORE::LINALG::SparseMatrix(*dofrowmap,numConstrID_,false,true));

  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("OffsetID", offset_id_);
  p.set("NumberofID", num_constr_id_);
  p.set("old disp", displast);
  p.set("new disp", disp);
  p.set("scaleStiffEntries", scStiff);
  p.set("scaleConstrMat", scConMat);
  p.set("vector curve factors", fact_);
  // Convert Epetra_Vector containing lagrange multipliers to an completely
  // redundant Epetra_vector since every element with the constraint condition needs them
  Teuchos::RCP<Epetra_Vector> lagrMultVecDense = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  CORE::LINALG::Export(*lagr_mult_vec_, *lagrMultVecDense);
  p.set("LagrMultVector", lagrMultVecDense);
  // Construct a redundant time curve factor and put it into parameter list
  Teuchos::RCP<Epetra_Vector> factredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  p.set("vector curve factors", factredundant);

  Teuchos::RCP<Epetra_Vector> actredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  Teuchos::RCP<Epetra_Vector> refbaseredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement", disp);
  volconstr3d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  areaconstr3d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  areaconstr2d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  volconstr3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  areaconstr3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);

  mpconplane3d_->SetConstrState("displacement", disp);
  mpconplane3d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  mpcnormcomp3d_->SetConstrState("displacement", disp);
  mpcnormcomp3d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  mpcnormcomp3dpen_->SetConstrState("displacement", disp);
  mpcnormcomp3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  mpconline2d_->SetConstrState("displacement", disp);
  mpconline2d_->Evaluate(p, stiff, constr_matrix_, fint, refbaseredundant, actredundant);
  // Export redundant vectors into distributed ones
  actvalues_->PutScalar(0.0);
  actvalues_->Export(*actredundant, *conimpo_, Add);
  Teuchos::RCP<Epetra_Vector> addrefbase = Teuchos::rcp(new Epetra_Vector(*constrmap_));
  addrefbase->Export(*refbaseredundant, *conimpo_, Add);
  refbasevalues_->Update(1.0, *addrefbase, 1.0);
  fact_->PutScalar(0.0);
  fact_->Export(*factredundant, *conimpo_, AbsMax);
  // ----------------------------------------------------
  // -----------include possible further constraints here
  // ----------------------------------------------------
  // Compute current reference volumes as elemetwise product of timecurvefactor and initialvalues
  referencevalues_->Multiply(1.0, *fact_, *refbasevalues_, 0.0);
  constrainterr_->Update(scConMat, *referencevalues_, -1.0 * scConMat, *actvalues_, 0.0);
  actdisc_->ClearState();
  // finalize the constraint matrix
  std::string label(constr_matrix_->Label());
  if (label == "CORE::LINALG::BlockSparseMatrixBase")
    constr_matrix_->Complete();
  else
    constr_matrix_->Complete(*constrmap_, *dofrowmap);
}

/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::ComputeError(double time, Teuchos::RCP<Epetra_Vector> disp)
{
  CheckIsInit();
  CheckIsSetup();

  std::vector<DRT::Condition*> constrcond(0);
  Teuchos::ParameterList p;
  p.set("total time", time);
  actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  CORE::LINALG::Export(*actvalues_, *actredundant);
  // Compute current values and assemble them to the completely redundant vector
  // We will always use the third systemvector for this purpose
  p.set("OffsetID", offset_id_);
  volconstr3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);
  areaconstr3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);
  areaconstr2d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);

  mpconplane3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);
  mpconplane3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);

  mpcnormcomp3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);
  mpcnormcomp3d_->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, actredundant);

  // Export redundant vectors into distributed ones
  actvalues_->PutScalar(0.0);
  actvalues_->Export(*actredundant, *conimpo_, Add);

  constrainterr_->Update(1.0, *referencevalues_, -1.0, *actvalues_, 0.0);
}


/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::ReadRestart(IO::DiscretizationReader& reader, const double& time)
{
  //  double uzawatemp = reader.ReadDouble("uzawaparameter");
  //  consolv_->SetUzawaParameter(uzawatemp);
  Teuchos::RCP<Epetra_Map> constrmap = GetConstraintMap();
  Teuchos::RCP<Epetra_Vector> tempvec = CORE::LINALG::CreateVector(*constrmap, true);
  reader.ReadVector(tempvec, "lagrmultiplier");
  SetLagrMultVector(tempvec);
  reader.ReadVector(tempvec, "refconval");
  SetRefBaseValues(tempvec, time);
}

/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::SetRefBaseValues(
    Teuchos::RCP<Epetra_Vector> newrefval, const double& time)
{
  volconstr3d_->Initialize(time);
  areaconstr3d_->Initialize(time);
  areaconstr2d_->Initialize(time);
  mpconplane3d_->Initialize(time);
  mpcnormcomp3d_->Initialize(time);
  mpconline2d_->Initialize(time);

  refbasevalues_->Update(1.0, *newrefval, 0.0);
}

/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::UpdateLagrMult(double factor)
{
  lagr_mult_vec_->Update(factor, *constrainterr_, 1.0);
  if (volconstr3d_->HaveConstraint())
  {
    std::vector<int> volconID = volconstr3d_->GetActiveCondID();
    for (unsigned int i = 0; i < volconID.size(); i++)
    {
      if (constrmap_->LID(int(i - offset_id_)) != -1)
      {
        std::cout << "Multiplier for Volume Constraint: " << volconID.at(i) << ":  "
                  << (*lagr_mult_vec_)[constrmap_->LID(int(i - offset_id_))] << '\n';
      }
    }
  }
}

void CONSTRAINTS::ConstrManager::Update() { lagr_mult_vec_old_->Update(1.0, *lagr_mult_vec_, 0.0); }

/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::UpdateLagrMult(Teuchos::RCP<Epetra_Vector> vect)
{
  lagr_mult_vec_->Update(1.0, *vect, 1.0);
}

void CONSTRAINTS::ConstrManager::UpdateTotLagrMult(Teuchos::RCP<Epetra_Vector> vect)
{
  lagr_mult_vec_->Update(1.0, *vect, 1.0, *lagr_mult_vec_old_, 0.0);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::ComputeMonitorValues(Teuchos::RCP<Epetra_Vector> disp)
{
  std::vector<DRT::Condition*> monitcond(0);
  monitorvalues_->PutScalar(0.0);
  Teuchos::ParameterList p;
  actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  p.set("OffsetID", min_monitor_id_);

  volmonitor3d_->Evaluate(p, actmonredundant);
  areamonitor3d_->Evaluate(p, actmonredundant);
  areamonitor2d_->Evaluate(p, actmonredundant);

  Epetra_Import monimpo(*monitormap_, *redmonmap_);
  monitorvalues_->Export(*actmonredundant, *monimpo_, Add);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::ComputeMonitorValues(Teuchos::RCP<const Epetra_Vector> disp)
{
  std::vector<DRT::Condition*> monitcond(0);
  monitorvalues_->PutScalar(0.0);
  Teuchos::ParameterList p;
  if (not actdisc_->DofRowMap()->SameAs(disp->Map()))
  {
    // build merged dof row map
    Teuchos::RCP<Epetra_Map> largemap =
        CORE::LINALG::MergeMap(*actdisc_->DofRowMap(), *constrmap_, false);

    CORE::LINALG::MapExtractor conmerger;
    conmerger.Setup(*largemap, Teuchos::rcp(actdisc_->DofRowMap(), false), constrmap_);
    actdisc_->SetState("displacement", conmerger.ExtractCondVector(disp));
  }
  else
    actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  p.set("OffsetID", min_monitor_id_);

  volmonitor3d_->Evaluate(p, actmonredundant);
  areamonitor3d_->Evaluate(p, actmonredundant);
  areamonitor2d_->Evaluate(p, actmonredundant);

  Epetra_Import monimpo(*monitormap_, *redmonmap_);
  monitorvalues_->Export(*actmonredundant, *monimpo_, Add);
}

/*----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::PrintMonitorValues() const
{
  if (num_monitor_id_ == 1)
    printf("Monitor value:\n");
  else if (num_monitor_id_ > 1)
    printf("Monitor values:\n");

  for (int i = 0; i < num_monitor_id_; ++i)
  {
    if ((*monitortypes_)[i] == 1.0)
    {
      printf("%2d (volume): %10.5e (%5.2f%% of initial value)\n", i + min_monitor_id_,
          abs((*monitorvalues_)[i]), ((*monitorvalues_)[i]) * 100 / ((*initialmonvalues_)[i]));
    }
    else if ((*monitortypes_)[i] == 2.0)
    {
      printf("%2d   (area): %10.5e (%5.2f%% of initial value)\n", i + min_monitor_id_,
          abs((*monitorvalues_)[i]), ((*monitorvalues_)[i]) * 100 / ((*initialmonvalues_)[i]));
    }
  }
}

void CONSTRAINTS::ConstrManager::BuildMoniType()
{
  Teuchos::ParameterList p1;
  // build distributed and redundant dummy monitor vector
  Teuchos::RCP<Epetra_Vector> dummymonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  Teuchos::RCP<Epetra_Vector> dummymondist = Teuchos::rcp(new Epetra_Vector(*monitormap_));
  p1.set("OffsetID", min_monitor_id_);

  // do the volumes
  volmonitor3d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  CORE::LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 1.0;
  }

  // do the area in 3D
  dummymonredundant->PutScalar(0.0);
  dummymondist->PutScalar(0.0);
  areamonitor3d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  CORE::LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 2.0;
  }

  // do the area in 2D
  dummymonredundant->PutScalar(0.0);
  dummymondist->PutScalar(0.0);
  areamonitor2d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  CORE::LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 3.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONSTRAINTS::ConstrManager::UseBlockMatrix(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps)
{
  // (re)allocate system matrix
  constr_matrix_ =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *domainmaps, *rangemaps, 81, false, true));
}

FOUR_C_NAMESPACE_CLOSE
