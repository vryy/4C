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
UTILS::ConstrManager::ConstrManager(RCP<DRT::Discretization> discr,
        RCP<Epetra_Vector> disp,
        ParameterList params):
actdisc_(discr)
{
  //Check, what kind of constraining boundary conditions there are
  numConstrID_=0;
  // Keep ParameterList p alive during initialization, so global information
  // over IDs as well as element results stored here can be used after all
  // constraints are evaluated
  ParameterList p;
  actdisc_->SetState("displacement",disp);
  minConstrID_=10000;
  int maxConstrID=0;
  volconstr3d_=rcp(new Constraint(actdisc_,"VolumeConstraint_3D",minConstrID_,maxConstrID));
  areaconstr3d_=rcp(new Constraint(actdisc_,"AreaConstraint_3D",minConstrID_,maxConstrID));
  areaconstr2d_=rcp(new Constraint(actdisc_,"AreaConstraint_2D",minConstrID_,maxConstrID));
  
  mpconplane3d_=rcp(new MPConstraint(actdisc_,"MPC_NodeOnPlane_3D",minConstrID_,maxConstrID));
  mpconline2d_=rcp(new MPConstraint(actdisc_,"MPC_NodeOnLine_2D",minConstrID_,maxConstrID));

  numConstrID_=max(maxConstrID-minConstrID_+1,0);
  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  haveconstraint_= (areaconstr3d_->HaveConstraint())||(volconstr3d_->HaveConstraint())||
        (areaconstr2d_->HaveConstraint())||(mpconplane3d_->HaveConstraint())||(mpconline2d_->HaveConstraint());
  if (haveconstraint_)
  {
    uzawaparam_=params.get<double>("uzawa parameter",1);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    //ManageIDs(p,minConstrID_,maxConstrID_,numConstrID_,MPCcondIDs);
    //initialize constrMatrix
    constrMatrix_=rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,true,true));
    //build Epetra_Map used as domainmap for constrMatrix and rowmap for result vectors 
    constrmap_=rcp(new Epetra_Map(numConstrID_,0,actdisc_->Comm()));
    //build an all reduced version of the constraintmap, since sometimes all processers
    //have to know all values of the constraints and Lagrange multipliers
    redconstrmap_ = LINALG::AllreduceEMap(*constrmap_);
    // sum up initial values
    initialvalues_=rcp(new Epetra_Vector(*constrmap_));
    RCP<Epetra_Vector> initialredundant = rcp(new Epetra_Vector(*redconstrmap_));
    LINALG::Export(*initialvalues_,*initialredundant);
    //Compute initial values and assemble them to the completely redundant vector
    //We will always use the third systemvector for this purpose
    p.set("MinID",minConstrID_);
    volconstr3d_->Evaluate(p,null,null,null,null,initialredundant);
    areaconstr3d_->Evaluate(p,null,null,null,null,initialredundant);
    areaconstr2d_->Evaluate(p,null,null,null,null,initialredundant);
    
    ImportResults(initialvalues_,initialredundant);
    
    mpconplane3d_->Evaluate(p,null,null,initialvalues_,null,null,true);
    mpconplane3d_->Evaluate(p,null,null,initialvalues_,null,null,true);

    //Initialize Lagrange Multiplicators, reference values and errors
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
  ParameterList p1;
  actdisc_->SetState("displacement",disp);
  minMonitorID_=10000;
  int maxMonitorID=0;
  volmonitor3d_=rcp(new Constraint(actdisc_,"VolumeMonitor_3D",minMonitorID_,maxMonitorID));
  areamonitor3d_=rcp(new Constraint(actdisc_,"AreaMonitor_3D",minMonitorID_,maxMonitorID));
  areamonitor2d_=rcp(new Constraint(actdisc_,"AreaMonitor_2D",minMonitorID_,maxMonitorID));
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  numMonitorID_=max(maxMonitorID-minMonitorID_+1,0);
  havemonitor_= (areamonitor3d_->HaveConstraint())||(volmonitor3d_->HaveConstraint())||(areamonitor2d_->HaveConstraint());
  if (havemonitor_)
  {
    //monitor values are only stored on processor zero since they are needed for output
    int nummyele=0;
    if (!actdisc_->Comm().MyPID())
    {
      nummyele=numMonitorID_;
    }
    monitormap_=rcp(new Epetra_Map(numMonitorID_,nummyele,0,actdisc_->Comm()));
    redmonmap_ = LINALG::AllreduceEMap(*monitormap_);
    monitorvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_->Scale(0.0);
    
    RCP<Epetra_Vector> initialmonredundant = rcp(new Epetra_Vector(*redmonmap_));
    LINALG::Export(*initialmonvalues_,*initialmonredundant);
    p1.set("MinID",minMonitorID_);
    
    volmonitor3d_->Evaluate(p1,null,null,null,null,initialmonredundant);
    areamonitor3d_->Evaluate(p1,null,null,null,null,initialmonredundant);
    areamonitor2d_->Evaluate(p1,null,null,null,null,initialmonredundant);

    ImportResults(initialmonvalues_,initialmonredundant);
    
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
        RCP<Epetra_Vector> disp,
        RCP<Epetra_Vector> fint,
        RCP<LINALG::SparseMatrix> stiff)
{
  // create the parameters for the discretization
  ParameterList p;
  vector<DRT::Condition*> constrcond(0);
  actvalues_->Scale(0.0);
  constrMatrix_->PutScalar(0.0);
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
  
  // other parameters that might be needed by the elements
  p.set("total time",time);
  p.set("MinID",minConstrID_);
  p.set("NumberofID",numConstrID_);
  // Convert Epetra_Vector constaining lagrange multipliers to an completely
  // redundant Epetra_vector since every element with the constraint condition needs them
  RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
  p.set("LagrMultVector",lagrMultVecDense);
  
  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);
 
  RCP<Epetra_Vector> actredundant = rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*actvalues_,*actredundant);

  volconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,null,actredundant);
  areaconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,null,actredundant);
  areaconstr2d_->Evaluate(p,stiff,constrMatrix_,fint,null,actredundant);
  
  mpconplane3d_->SetConstrState("displacement",disp);
  mpconplane3d_->Evaluate(p,stiff,constrMatrix_,fint,null,actredundant);
  mpconline2d_->SetConstrState("displacement",disp);
  mpconline2d_->Evaluate(p,stiff,constrMatrix_,fint,null,actredundant);
  
  ImportResults(actvalues_,actredundant);
  
  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  SynchronizeMinConstraint(p,fact_,"LoadCurveFactor");
  //Compute current referencevolumes as elemetwise product of timecurvefactor and initialvalues
  referencevalues_->Multiply(1.0,*fact_,*initialvalues_,0.0);
  constrainterr_->Update(1.0,*referencevalues_,-1.0,*actvalues_,0.0);
  actdisc_->ClearState();
  // finalize the constraint matrix
  constrMatrix_->Complete(*constrmap_,*dofrowmap);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Compute difference between current and prescribed values.              |
*-----------------------------------------------------------------------*/
void UTILS::ConstrManager::ComputeError(double time,RCP<Epetra_Vector> disp)
{
    vector<DRT::Condition*> constrcond(0);
    actvalues_->Scale(0.0);
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
    
    ImportResults(actvalues_,actredundant);

    constrainterr_->Update(1.0,*referencevalues_,-1.0,*actvalues_,0.0);
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
  LINALG::Export(*monitorvalues_,*actmonredundant);
  p.set("MinID",minMonitorID_);
  
  volmonitor3d_->Evaluate(p,null,null,null,null,actmonredundant);
  areamonitor3d_->Evaluate(p,null,null,null,null,actmonredundant);
  areamonitor2d_->Evaluate(p,null,null,null,null,actmonredundant);
  
  ImportResults(monitorvalues_,actmonredundant);
  
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

/*----------------------------------------------------------------------*
 |(private)                                                 tk 07/08    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraints by summing them up                                        |
 *----------------------------------------------------------------------*/
void UTILS::ConstrManager::ImportResults
(
  RCP<Epetra_Vector>& vect_dist,
  RCP<Epetra_Vector>& vect_redu
)
{
  vector<int> gids;
  for (int i = 0; i < (vect_redu->MyLength()); ++i)
  {
    gids.push_back(i);
  }
  vector<double> currval(vect_redu->MyLength(),0);
  actdisc_->Comm().SumAll(&((*vect_redu)[0]),&(currval[0]),vect_redu->MyLength());
  vect_dist->SumIntoGlobalValues(vect_redu->MyLength(),&(currval[0]),&(gids[0]));

  return;
}


/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraint by finding minimum value                                   |
 *----------------------------------------------------------------------*/
void UTILS::ConstrManager::SynchronizeMinConstraint(ParameterList& params,
                                RCP<Epetra_Vector>& vect,
                                const char* resultstring)
{

  for (int i = 0; i < numConstrID_; ++i)
  {
    char valname[30];
    sprintf(valname,"%s %d",resultstring,i+minConstrID_);
    double currval=1;
    actdisc_->Comm().MinAll(&(params.get(valname,1.0)),&(currval),1);
    vect->ReplaceGlobalValues(1,&currval,&i);
  }

  return;
}

#endif
