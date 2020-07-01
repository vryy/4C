/*----------------------------------------------------------------------*/
/*! \file

\brief Class controlling constraints and containing the necessary data, code originally by Thomas
Kloeppel


\level 2


*----------------------------------------------------------------------*/


#include <iostream>

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_io/io.H"

#include "constraint_manager.H"
#include "constraint.H"
#include "constraintpenalty.H"
#include "monitor.H"
#include "multipointconstraint2.H"
#include "multipointconstraint3.H"
#include "multipointconstraint3penalty.H"
#include "constraintdofset.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
UTILS::ConstrManager::ConstrManager()
    : offsetID_(-1),
      maxConstrID_(0),
      numConstrID_(-1),
      numMonitorID_(-1),
      minMonitorID_(-1),
      haveconstraint_(false),
      havelagrconstr_(false),
      havepenaconstr_(false),
      havemonitor_(false),
      uzawaparam_(0.0),
      issetup_(false),
      isinit_(false)

{
  // Keep constructor empty !
  return;
}


/*----------------------------------------------------------------------*
 |  initialize this class                                   rauch 09/16 |
 *----------------------------------------------------------------------*/
void UTILS::ConstrManager::Init(
    Teuchos::RCP<DRT::Discretization> discr, Teuchos::ParameterList params)
{
  SetIsSetup(false);

  // set pointer to discretization
  actdisc_ = discr;

  //----------------------------------------------------------------------------
  //---------------------------------------------------------Constraint Conditions!

  // constructors of constraints increment number of constraints defined and the minimum
  // ConditionID read so far.
  numConstrID_ = 0;
  offsetID_ = 10000;
  // Check, what kind of constraining boundary conditions there are
  volconstr3d_ =
      Teuchos::rcp(new Constraint(actdisc_, "VolumeConstraint_3D", offsetID_, maxConstrID_));
  areaconstr3d_ =
      Teuchos::rcp(new Constraint(actdisc_, "AreaConstraint_3D", offsetID_, maxConstrID_));
  areaconstr2d_ =
      Teuchos::rcp(new Constraint(actdisc_, "AreaConstraint_2D", offsetID_, maxConstrID_));
  mpconline2d_ =
      Teuchos::rcp(new MPConstraint2(actdisc_, "MPC_NodeOnLine_2D", offsetID_, maxConstrID_));
  mpconplane3d_ =
      Teuchos::rcp(new MPConstraint3(actdisc_, "MPC_NodeOnPlane_3D", offsetID_, maxConstrID_));
  mpcnormcomp3d_ =
      Teuchos::rcp(new MPConstraint3(actdisc_, "MPC_NormalComponent_3D", offsetID_, maxConstrID_));

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
  return;
}


/*----------------------------------------------------------------------*
 |  setup this class                                        rauch 09/16 |
 *----------------------------------------------------------------------*/
void UTILS::ConstrManager::Setup(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::ParameterList params)
{
  CheckIsInit();

  if (haveconstraint_)
  {
    numConstrID_ = std::max(maxConstrID_ - offsetID_ + 1, 0);
    constrdofset_ = Teuchos::rcp(new ConstraintDofSet());
    constrdofset_->AssignDegreesOfFreedom(actdisc_, numConstrID_, 0);
    offsetID_ -= constrdofset_->FirstGID();
    Teuchos::ParameterList p;
    uzawaparam_ = params.get<double>("uzawa parameter", 1);
    double time = params.get<double>("total time", 0.0);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    // initialize constrMatrix
    constrMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, numConstrID_, false, true));
    // build Epetra_Map used as domainmap for constrMatrix and rowmap for result vectors
    constrmap_ = Teuchos::rcp(new Epetra_Map(*(constrdofset_->DofRowMap())));
    // build an all reduced version of the constraintmap, since sometimes all processors
    // have to know all values of the constraints and Lagrange multipliers
    redconstrmap_ = LINALG::AllreduceEMap(*constrmap_);
    // importer
    conimpo_ = Teuchos::rcp(new Epetra_Export(*redconstrmap_, *constrmap_));
    // sum up initial values
    refbasevalues_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
    Teuchos::RCP<Epetra_Vector> refbaseredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
    // Compute initial values and assemble them to the completely redundant vector
    // We will always use the third systemvector for this purpose
    p.set("OffsetID", offsetID_);
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
    lagrMultVec_ = Teuchos::rcp(new Epetra_Vector(*constrmap_, true));
    lagrMultVecOld_ = Teuchos::rcp(new Epetra_Vector(*constrmap_, true));
    fact_ = Teuchos::rcp(new Epetra_Vector(*constrmap_));
  }
  //----------------------------------------------------------------------------
  //---------------------------------------------------------Monitor Conditions!
  actdisc_->SetState("displacement", disp);
  minMonitorID_ = 10000;
  int maxMonitorID = 0;
  volmonitor3d_ =
      Teuchos::rcp(new Monitor(actdisc_, "VolumeMonitor_3D", minMonitorID_, maxMonitorID));
  areamonitor3d_ =
      Teuchos::rcp(new Monitor(actdisc_, "AreaMonitor_3D", minMonitorID_, maxMonitorID));
  areamonitor2d_ =
      Teuchos::rcp(new Monitor(actdisc_, "AreaMonitor_2D", minMonitorID_, maxMonitorID));
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  numMonitorID_ = std::max(maxMonitorID - minMonitorID_ + 1, 0);
  havemonitor_ = (areamonitor3d_->HaveMonitor()) || (volmonitor3d_->HaveMonitor()) ||
                 (areamonitor2d_->HaveMonitor());
  if (havemonitor_)
  {
    Teuchos::ParameterList p1;
    // monitor values are only stored on processor zero since they are only needed for output
    int nummyele = 0;
    if (!actdisc_->Comm().MyPID())
    {
      nummyele = numMonitorID_;
    }
    // initialize maps and importer
    monitormap_ = Teuchos::rcp(new Epetra_Map(numMonitorID_, nummyele, 0, actdisc_->Comm()));
    redmonmap_ = LINALG::AllreduceEMap(*monitormap_);
    monimpo_ = Teuchos::rcp(new Epetra_Export(*redmonmap_, *monitormap_));
    monitorvalues_ = Teuchos::rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_ = Teuchos::rcp(new Epetra_Vector(*monitormap_));

    Teuchos::RCP<Epetra_Vector> initialmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
    p1.set("OffsetID", minMonitorID_);
    volmonitor3d_->Evaluate(p1, initialmonredundant);
    areamonitor3d_->Evaluate(p1, initialmonredundant);
    areamonitor2d_->Evaluate(p1, initialmonredundant);

    // Export redundant vector into distributed one
    initialmonvalues_->Export(*initialmonredundant, *monimpo_, Add);
    monitortypes_ = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
    BuildMoniType();
  }

  SetIsSetup(true);
  return;
}


/*----------------------------------------------------------------------*
|(public)                                                       tk 11/07|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::EvaluateForceStiff(const double time,
    Teuchos::RCP<const Epetra_Vector> displast, Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> fint, Teuchos::RCP<LINALG::SparseOperator> stiff,
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
  constrMatrix_
      ->Zero();  //=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,false,true));

  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("OffsetID", offsetID_);
  p.set("NumberofID", numConstrID_);
  p.set("old disp", displast);
  p.set("new disp", disp);
  p.set("scaleStiffEntries", scStiff);
  p.set("scaleConstrMat", scConMat);
  p.set("vector curve factors", fact_);
  // Convert Epetra_Vector containing lagrange multipliers to an completely
  // redundant Epetra_vector since every element with the constraint condition needs them
  Teuchos::RCP<Epetra_Vector> lagrMultVecDense = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*lagrMultVec_, *lagrMultVecDense);
  p.set("LagrMultVector", lagrMultVecDense);
  // Construct a redundant time curve factor and put it into parameter list
  Teuchos::RCP<Epetra_Vector> factredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  p.set("vector curve factors", factredundant);

  Teuchos::RCP<Epetra_Vector> actredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  Teuchos::RCP<Epetra_Vector> refbaseredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement", disp);
  volconstr3d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
  areaconstr3d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
  areaconstr2d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
  volconstr3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  areaconstr3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);

  mpconplane3d_->SetConstrState("displacement", disp);
  mpconplane3d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
  mpcnormcomp3d_->SetConstrState("displacement", disp);
  mpcnormcomp3d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
  mpcnormcomp3dpen_->SetConstrState("displacement", disp);
  mpcnormcomp3dpen_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  mpconline2d_->SetConstrState("displacement", disp);
  mpconline2d_->Evaluate(p, stiff, constrMatrix_, fint, refbaseredundant, actredundant);
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
  std::string label(constrMatrix_->Label());
  if (label == "LINALG::BlockSparseMatrixBase")
    constrMatrix_->Complete();
  else
    constrMatrix_->Complete(*constrmap_, *dofrowmap);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Compute difference between current and prescribed values.              |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeError(double time, Teuchos::RCP<Epetra_Vector> disp)
{
  CheckIsInit();
  CheckIsSetup();

  std::vector<DRT::Condition*> constrcond(0);
  Teuchos::ParameterList p;
  p.set("total time", time);
  actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actredundant = Teuchos::rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*actvalues_, *actredundant);
  // Compute current values and assemble them to the completely redundant vector
  // We will always use the third systemvector for this purpose
  p.set("OffsetID", offsetID_);
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
  return;
}


/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ReadRestart(IO::DiscretizationReader& reader, const double& time)
{
  //  double uzawatemp = reader.ReadDouble("uzawaparameter");
  //  consolv_->SetUzawaParameter(uzawatemp);
  Teuchos::RCP<Epetra_Map> constrmap = GetConstraintMap();
  Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*constrmap, true);
  reader.ReadVector(tempvec, "lagrmultiplier");
  SetLagrMultVector(tempvec);
  reader.ReadVector(tempvec, "refconval");
  SetRefBaseValues(tempvec, time);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 08/08|
|Reset reference base values for restart                                |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::SetRefBaseValues(
    Teuchos::RCP<Epetra_Vector> newrefval, const double& time)
{
  volconstr3d_->Initialize(time);
  areaconstr3d_->Initialize(time);
  areaconstr2d_->Initialize(time);
  mpconplane3d_->Initialize(time);
  mpcnormcomp3d_->Initialize(time);
  mpconline2d_->Initialize(time);

  refbasevalues_->Update(1.0, *newrefval, 0.0);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add scaled error of constraint to Lagrange multiplier.                 |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::UpdateLagrMult(double factor)
{
  lagrMultVec_->Update(factor, *constrainterr_, 1.0);
  if (volconstr3d_->HaveConstraint())
  {
    std::vector<int> volconID = volconstr3d_->GetActiveCondID();
    for (unsigned int i = 0; i < volconID.size(); i++)
    {
      if (constrmap_->LID(int(i - offsetID_)) != -1)
      {
        std::cout << "Multiplier for Volume Constraint: " << volconID.at(i) << ":  "
                  << (*lagrMultVec_)[constrmap_->LID(int(i - offsetID_))] << std::endl;
      }
    }
  }

  return;
}

void UTILS::ConstrManager::Update() { lagrMultVecOld_->Update(1.0, *lagrMultVec_, 0.0); }

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add Lagrange increment to Lagrange multiplier.                         |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::UpdateLagrMult(Teuchos::RCP<Epetra_Vector> vect)
{
  lagrMultVec_->Update(1.0, *vect, 1.0);
  return;
}

void UTILS::ConstrManager::UpdateTotLagrMult(Teuchos::RCP<Epetra_Vector> vect)
{
  lagrMultVec_->Update(1.0, *vect, 1.0, *lagrMultVecOld_, 0.0);
  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeMonitorValues(Teuchos::RCP<Epetra_Vector> disp)
{
  std::vector<DRT::Condition*> monitcond(0);
  monitorvalues_->PutScalar(0.0);
  Teuchos::ParameterList p;
  actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  p.set("OffsetID", minMonitorID_);

  volmonitor3d_->Evaluate(p, actmonredundant);
  areamonitor3d_->Evaluate(p, actmonredundant);
  areamonitor2d_->Evaluate(p, actmonredundant);

  Epetra_Import monimpo(*monitormap_, *redmonmap_);
  monitorvalues_->Export(*actmonredundant, *monimpo_, Add);

  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeMonitorValues(Teuchos::RCP<const Epetra_Vector> disp)
{
  std::vector<DRT::Condition*> monitcond(0);
  monitorvalues_->PutScalar(0.0);
  Teuchos::ParameterList p;
  if (not actdisc_->DofRowMap()->SameAs(disp->Map()))
  {
    // build merged dof row map
    Teuchos::RCP<Epetra_Map> largemap =
        LINALG::MergeMap(*actdisc_->DofRowMap(), *constrmap_, false);

    LINALG::MapExtractor conmerger;
    conmerger.Setup(*largemap, Teuchos::rcp(actdisc_->DofRowMap(), false), constrmap_);
    actdisc_->SetState("displacement", conmerger.ExtractCondVector(disp));
  }
  else
    actdisc_->SetState("displacement", disp);

  Teuchos::RCP<Epetra_Vector> actmonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  p.set("OffsetID", minMonitorID_);

  volmonitor3d_->Evaluate(p, actmonredundant);
  areamonitor3d_->Evaluate(p, actmonredundant);
  areamonitor2d_->Evaluate(p, actmonredundant);

  Epetra_Import monimpo(*monitormap_, *redmonmap_);
  monitorvalues_->Export(*actmonredundant, *monimpo_, Add);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Print monitored values                                                 |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::PrintMonitorValues() const
{
  if (numMonitorID_ == 1)
    printf("Monitor value:\n");
  else if (numMonitorID_ > 1)
    printf("Monitor values:\n");

  for (int i = 0; i < numMonitorID_; ++i)
  {
    if ((*monitortypes_)[i] == 1.0)
      printf("%2d (volume): %10.5e (%5.2f%% of initial value)\n", i + minMonitorID_,
          abs((*monitorvalues_)[i]), ((*monitorvalues_)[i]) * 100 / ((*initialmonvalues_)[i]));
    else if ((*monitortypes_)[i] == 2.0)
      printf("%2d   (area): %10.5e (%5.2f%% of initial value)\n", i + minMonitorID_,
          abs((*monitorvalues_)[i]), ((*monitorvalues_)[i]) * 100 / ((*initialmonvalues_)[i]));
    else if ((*monitortypes_)[i] == 2.0)
      printf("%2d   (area): %10.5e (%5.2f%% of initial value)\n", i + minMonitorID_,
          abs((*monitorvalues_)[i]), ((*monitorvalues_)[i]) * 100 / ((*initialmonvalues_)[i]));
  }

  return;
}

void UTILS::ConstrManager::BuildMoniType()
{
  Teuchos::ParameterList p1;
  // build distributed and redundant dummy monitor vector
  Teuchos::RCP<Epetra_Vector> dummymonredundant = Teuchos::rcp(new Epetra_Vector(*redmonmap_));
  Teuchos::RCP<Epetra_Vector> dummymondist = Teuchos::rcp(new Epetra_Vector(*monitormap_));
  p1.set("OffsetID", minMonitorID_);

  // do the volumes
  volmonitor3d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 1.0;
  }

  // do the area in 3D
  dummymonredundant->Scale(0.0);
  dummymondist->Scale(0.0);
  areamonitor3d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 2.0;
  }

  // do the area in 2D
  dummymonredundant->Scale(0.0);
  dummymondist->Scale(0.0);
  areamonitor2d_->Evaluate(p1, dummymonredundant);
  // Export redundant vector into distributed one
  dummymondist->Export(*dummymonredundant, *monimpo_, Add);
  // Now export back
  LINALG::Export(*dummymondist, *dummymonredundant);
  for (int i = 0; i < dummymonredundant->MyLength(); i++)
  {
    if ((*dummymonredundant)[i] != 0.0) (*monitortypes_)[i] = 3.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void UTILS::ConstrManager::UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
{
  // (re)allocate system matrix
  constrMatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *domainmaps, *rangemaps, 81, false, true));

  return;
}
