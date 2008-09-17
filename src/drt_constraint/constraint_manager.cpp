/*!----------------------------------------------------------------------
\file constraint_manager.cpp

\brief Class controlling constraints and containing the necessary data
<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraint_manager.H"
#include "constraint_element.H"
#include "mpcdofset.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
UTILS::ConstrManager::ConstrManager
(
  RCP<DRT::Discretization> discr,
  RCP<Epetra_Vector> disp,
  ParameterList params):
actdisc_(discr)
{
  //----------------------------------------------------------------------------
  //---------------------------------------------------------Constraint Conditions!
  
  //Check, what kind of constraining boundary conditions there are
  numConstrID_=0;
  // Keep ParameterList p alive during initialization, so global information
  // over IDs as well as element results stored here can be used after all
  // constraints are evaluated
  const Epetra_Map* dofrowmap1 = actdisc_->DofRowMap();
  int offset;
  offset=dofrowmap1->MaxAllGID()+1;
  minConstrID_=10000;
  int maxConstrID=0;
  volconstr3d_=rcp(new Constraint(actdisc_,"VolumeConstraint_3D",minConstrID_,maxConstrID));
  areaconstr3d_=rcp(new Constraint(actdisc_,"AreaConstraint_3D",minConstrID_,maxConstrID));
  areaconstr2d_=rcp(new Constraint(actdisc_,"AreaConstraint_2D",minConstrID_,maxConstrID));
  mpconline2d_=rcp(new MPConstraint2(actdisc_,"MPC_NodeOnLine_2D",minConstrID_,maxConstrID));
  mpconplane3d_=rcp(new MPConstraint3(actdisc_,"MPC_NodeOnPlane_3D",minConstrID_,maxConstrID));

  numConstrID_=max(maxConstrID-minConstrID_+1,0);
  minConstrID_-=offset;
  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  haveconstraint_= (areaconstr3d_->HaveConstraint())||(volconstr3d_->HaveConstraint())||
        (areaconstr2d_->HaveConstraint())||(mpconplane3d_->HaveConstraint())||(mpconline2d_->HaveConstraint());
  if (haveconstraint_)
  {
    ParameterList p;
    uzawaparam_ = params.get<double>("uzawa parameter",1);
    double time = params.get<double>("total time"      ,0.0);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    //initialize constrMatrix
    constrMatrix_=rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,true,true));
    //build Epetra_Map used as domainmap for constrMatrix and rowmap for result vectors 
    constrmap_=rcp(new Epetra_Map(numConstrID_,offset,actdisc_->Comm()));
    //build an all reduced version of the constraintmap, since sometimes all processors
    //have to know all values of the constraints and Lagrange multipliers
    redconstrmap_ = LINALG::AllreduceEMap(*constrmap_);
    // importer
    conimpo_ = rcp (new Epetra_Export(*redconstrmap_,*constrmap_));
    // sum up initial values
    refbasevalues_=rcp(new Epetra_Vector(*constrmap_));
    RCP<Epetra_Vector> refbaseredundant = rcp(new Epetra_Vector(*redconstrmap_));
    //Compute initial values and assemble them to the completely redundant vector
    //We will always use the third systemvector for this purpose
    p.set("MinID",minConstrID_);
    p.set("total time",time);
    actdisc_->SetState("displacement",disp);
    volconstr3d_->Initialize(p,refbaseredundant);
    areaconstr3d_->Initialize(p,refbaseredundant);
    areaconstr2d_->Initialize(p,refbaseredundant);
    
    mpconline2d_->SetConstrState("displacement",disp);
    mpconline2d_->Initialize(p,refbasevalues_);
    mpconplane3d_->SetConstrState("displacement",disp);
    mpconplane3d_->Initialize(p,refbaseredundant);
    
    // Export redundant vector into distributed one
    refbasevalues_ -> Export(*refbaseredundant,*conimpo_,Add);
    
    //Initialize Lagrange Multipliers, reference values and errors
    actdisc_->ClearState();
    referencevalues_=rcp(new Epetra_Vector(*constrmap_));
    actvalues_=rcp(new Epetra_Vector(*constrmap_));
    actvalues_->Scale(0.0);
    constrainterr_=rcp(new Epetra_Vector(*constrmap_));
    lagrMultVec_=rcp(new Epetra_Vector(*constrmap_));
    lagrMultVec_->Scale(0.0);
    fact_=rcp(new Epetra_Vector(*constrmap_));
  }
  //----------------------------------------------------------------------------
  //---------------------------------------------------------Monitor Conditions!
  actdisc_->SetState("displacement",disp);
  minMonitorID_=10000;
  int maxMonitorID=0;
  volmonitor3d_=rcp(new Monitor(actdisc_,"VolumeMonitor_3D",minMonitorID_,maxMonitorID));
  areamonitor3d_=rcp(new Monitor(actdisc_,"AreaMonitor_3D",minMonitorID_,maxMonitorID));
  areamonitor2d_=rcp(new Monitor(actdisc_,"AreaMonitor_2D",minMonitorID_,maxMonitorID));
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  numMonitorID_=max(maxMonitorID-minMonitorID_+1,0);
  havemonitor_= (areamonitor3d_->HaveMonitor())||(volmonitor3d_->HaveMonitor())||(areamonitor2d_->HaveMonitor());
  if (havemonitor_)
  {

    ParameterList p1;
    //monitor values are only stored on processor zero since they are only needed for output
    int nummyele=0;
    if (!actdisc_->Comm().MyPID())
    {
      nummyele=numMonitorID_;
    }
    // initialize maps and importer
    monitormap_=rcp(new Epetra_Map(numMonitorID_,nummyele,0,actdisc_->Comm()));
    redmonmap_ = LINALG::AllreduceEMap(*monitormap_);
    monimpo_ = rcp (new Epetra_Export(*redmonmap_,*monitormap_));
    monitorvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_=rcp(new Epetra_Vector(*monitormap_));
    
    RCP<Epetra_Vector> initialmonredundant = rcp(new Epetra_Vector(*redmonmap_));
    LINALG::Export(*initialmonvalues_,*initialmonredundant);
    p1.set("MinID",minMonitorID_);
    volmonitor3d_->Evaluate(p1,initialmonredundant);
    areamonitor3d_->Evaluate(p1,initialmonredundant);
    areamonitor2d_->Evaluate(p1,initialmonredundant);

    // Export redundant vector into distributed one
    initialmonvalues_->Export(*initialmonredundant,*monimpo_,Add);
  }
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 11/07|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::StiffnessAndInternalForces(
        const double time,
        RCP<Epetra_Vector> displast,
        RCP<Epetra_Vector> disp,
        RCP<Epetra_Vector> fint,
        RCP<LINALG::SparseMatrix> stiff,
        ParameterList scalelist)
{
  double scStiff = scalelist.get("scaleStiffEntries",1.0);
  double scConMat = scalelist.get("scaleConstrMat",1.0);
  
  // create the parameters for the discretization
  ParameterList p;
  vector<DRT::Condition*> constrcond(0);
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
  constrMatrix_->Scale(0.0);//rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,true,true));
    
  // other parameters that might be needed by the elements
  p.set("total time",time);
  p.set("MinID",minConstrID_);
  p.set("NumberofID",numConstrID_);
  p.set("old disp",displast);
  p.set("new disp",disp);
  p.set("scaleStiffEntries",scStiff);
  p.set("scaleConstrMat",scConMat);
  p.set("vector curve factors",fact_);
  // Convert Epetra_Vector containing lagrange multipliers to an completely
  // redundant Epetra_vector since every element with the constraint condition needs them
  RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
  p.set("LagrMultVector",lagrMultVecDense);
  // Construct a redundant time curve factor and put it into parameter list
  RCP<Epetra_Vector> factredundant = rcp(new Epetra_Vector(*redconstrmap_));
  p.set("vector curve factors",factredundant);
  
  RCP<Epetra_Vector> actredundant = rcp(new Epetra_Vector(*redconstrmap_));
  RCP<Epetra_Vector> refbaseredundant = rcp(new Epetra_Vector(*redconstrmap_));
  
  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);
  volconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,refbaseredundant,actredundant);
  areaconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,refbaseredundant,actredundant);
  areaconstr2d_->Evaluate(p,stiff,constrMatrix_,fint,refbaseredundant,actredundant);
  
  mpconplane3d_->SetConstrState("displacement",disp);
  mpconplane3d_->Evaluate(p,stiff,constrMatrix_,fint,refbaseredundant,actredundant);
  mpconline2d_->SetConstrState("displacement",disp);
  mpconline2d_->Evaluate(p,stiff,constrMatrix_,fint,refbaseredundant,actredundant);
  
  // Export redundant vectors into distributed ones
  actvalues_->Scale(0.0);
  actvalues_->Export(*actredundant,*conimpo_,Add);
  RCP<Epetra_Vector> addrefbase = rcp(new Epetra_Vector(*constrmap_));
  addrefbase->Export(*refbaseredundant,*conimpo_,Add);
  refbasevalues_->Update(1.0,*addrefbase,1.0);
  fact_->Scale(0.0);
  fact_->Export(*factredundant,*conimpo_,Insert);
  // ----------------------------------------------------
  // -----------include possible further constraints here
  // ----------------------------------------------------
  // Compute current reference volumes as elemetwise product of timecurvefactor and initialvalues
  referencevalues_->Multiply(1.0,*fact_,*refbasevalues_,0.0);  
  constrainterr_->Update(scConMat,*referencevalues_,-1.0*scConMat,*actvalues_,0.0);
  actdisc_->ClearState();
  // finalize the constraint matrix
  constrMatrix_->Complete(*constrmap_,*dofrowmap);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Compute difference between current and prescribed values.              |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeError(double time, RCP<Epetra_Vector> disp)
{
    vector<DRT::Condition*> constrcond(0);
    ParameterList p;
    p.set("total time",time);
    actdisc_->SetState("displacement",disp);
    
    RCP<Epetra_Vector> actredundant = rcp(new Epetra_Vector(*redconstrmap_));
    LINALG::Export(*actvalues_,*actredundant);
    //Compute current values and assemble them to the completely redundant vector
    //We will always use the third systemvector for this purpose
    p.set("MinID",minConstrID_);
    volconstr3d_->Evaluate(p,null,null,null,null,actredundant);
    areaconstr3d_->Evaluate(p,null,null,null,null,actredundant);
    areaconstr2d_->Evaluate(p,null,null,null,null,actredundant);
     
    mpconplane3d_->Evaluate(p,null,null,null,null,actredundant);
    mpconplane3d_->Evaluate(p,null,null,null,null,actredundant);
    
    // Export redundant vectors into distributed ones
    actvalues_->Scale(0.0);
    actvalues_->Export(*actredundant,*conimpo_,Add);

    constrainterr_->Update(1.0,*referencevalues_,-1.0,*actvalues_,0.0);
    return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 08/08|
|Reset reference base values for restart                                |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::SetRefBaseValues(RCP<Epetra_Vector> newrefval,const double& time)
{
  volconstr3d_->Initialize(time);
  areaconstr3d_->Initialize(time);
  areaconstr2d_->Initialize(time);
  mpconplane3d_->Initialize(time);
  mpconline2d_->Initialize(time);

  refbasevalues_->Update(1.0, *newrefval,0.0);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add scaled error of constraint to Lagrange multiplier.                 |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::UpdateLagrMult(double factor)
{
  lagrMultVec_->Update(factor,*constrainterr_,1.0);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add Lagrange increment to Lagrange multiplier.                         |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::UpdateLagrMult(RCP<Epetra_Vector> vect)
{
  lagrMultVec_->Update(1.0,*vect,1.0);
  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeMonitorValues(RCP<Epetra_Vector> disp)
{
  vector<DRT::Condition*> monitcond(0);
  monitorvalues_->Scale(0.0);
  ParameterList p;
  actdisc_->SetState("displacement",disp);
  
  RCP<Epetra_Vector> actmonredundant = rcp(new Epetra_Vector(*redmonmap_));
  p.set("MinID",minMonitorID_);
  
  volmonitor3d_->Evaluate(p,actmonredundant);
  areamonitor3d_->Evaluate(p,actmonredundant);
  areamonitor2d_->Evaluate(p,actmonredundant);
 
  Epetra_Import monimpo(*monitormap_,*redmonmap_);
  monitorvalues_->Export(*actmonredundant,*monimpo_,Add);
  
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Print monitored values                                                 |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::PrintMonitorValues()
{
  for (int i = 0; i < numMonitorID_; ++i)
  {
    printf("Monitor value %2d: %10.5e (%5.2f%% of initial value)\n",i+minMonitorID_,(*monitorvalues_)[i],((*monitorvalues_)[i])*100/((*initialmonvalues_)[i]));
  }

  return;
}

#endif
