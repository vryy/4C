/*!----------------------------------------------------------------------
\file drt_constraint_manager.cpp

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
ConstrManager::ConstrManager(RCP<DRT::Discretization> discr,
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
  maxConstrID_=0;
  volconstr3d_=rcp(new Constraint(actdisc_,"VolumeConstraint_3D",minConstrID_,maxConstrID_));
  areaconstr3d_=rcp(new Constraint(actdisc_,"AreaConstraint_3D",minConstrID_,maxConstrID_));
  areaconstr2d_=rcp(new Constraint(actdisc_,"AreaConstraint_2D",minConstrID_,maxConstrID_));
  
  mpconplane3d_=rcp(new MPConstraint(actdisc_,"MPC_NodeOnPlane_3D",minConstrID_,maxConstrID_));
  mpconline2d_=rcp(new MPConstraint(actdisc_,"MPC_NodeOnLine_2D",minConstrID_,maxConstrID_));

  numConstrID_=max(maxConstrID_-minConstrID_+1,0);
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
    //Set initial values to computed volumes and areas and to amplitudes of MPC
    volconstr3d_->Evaluate(p,null,null,null,null,null);
    areaconstr3d_->Evaluate(p,null,null,null,null,null);
    areaconstr2d_->Evaluate(p,null,null,null,null,null);
    SynchronizeSumConstraint(p,initialvalues_,"computed volume",numConstrID_,minConstrID_);
    SynchronizeSumConstraint(p,initialvalues_,"computed area",numConstrID_,minConstrID_);
    p.set("MinID",minConstrID_);
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
  maxMonitorID_=0;
  volmonitor3d_=rcp(new Constraint(actdisc_,"VolumeMonitor_3D",minMonitorID_,maxMonitorID_));
  areamonitor3d_=rcp(new Constraint(actdisc_,"AreaMonitor_3D",minMonitorID_,maxMonitorID_));
  areamonitor2d_=rcp(new Constraint(actdisc_,"AreaMonitor_2D",minMonitorID_,maxMonitorID_));
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  numMonitorID_=max(maxMonitorID_-minMonitorID_+1,0);
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
    monitorvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_->Scale(0.0);
    volmonitor3d_->Evaluate(p1,null,null,null,null,null);
    areamonitor3d_->Evaluate(p1,null,null,null,null,null);
    areamonitor2d_->Evaluate(p1,null,null,null,null,null);

    SynchronizeSumConstraint(p1,initialmonvalues_,"computed volume",numMonitorID_,minMonitorID_);
    SynchronizeSumConstraint(p1,initialmonvalues_,"computed area",numMonitorID_,minMonitorID_);
  }
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 11/07|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
*-----------------------------------------------------------------------*/
void ConstrManager::StiffnessAndInternalForces(
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
  p.set("MaxID",maxConstrID_);
  p.set("NumberofID",numConstrID_);
  // Convert Epetra_Vector constaining lagrange multipliers to an completely
  // redundant Epetra_vector since every element with the constraint condition needs them
  RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*redconstrmap_));
  LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
  p.set("LagrMultVector",lagrMultVecDense);
  
  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);
  
  volconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,null,null);
  SynchronizeSumConstraint(p,actvalues_,"computed volume",numConstrID_,minConstrID_);
  areaconstr3d_->Evaluate(p,stiff,constrMatrix_,fint,null,null);
  areaconstr2d_->Evaluate(p,stiff,constrMatrix_,fint,null,null);
  SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
  
  mpconplane3d_->SetConstrState("displacement",disp);
  mpconplane3d_->Evaluate(p,stiff,constrMatrix_,fint,null,null);
  mpconline2d_->SetConstrState("displacement",disp);
  mpconline2d_->Evaluate(p,stiff,constrMatrix_,fint,null,null);
  SynchronizeSumConstraint(p,actvalues_,"computed MPC value",numConstrID_,minConstrID_);
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
void ConstrManager::ComputeError(double time,RCP<Epetra_Vector> disp)
{
    vector<DRT::Condition*> constrcond(0);
    actvalues_->Scale(0.0);
    ParameterList p;
    p.set("total time",time);
    actdisc_->SetState("displacement",disp);
    
    volconstr3d_->Evaluate(p,null,null,null,null,null);
    SynchronizeSumConstraint(p,actvalues_,"computed volume",numConstrID_,minConstrID_);
    areaconstr3d_->Evaluate(p,null,null,null,null,null);
    areaconstr2d_->Evaluate(p,null,null,null,null,null);
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);

    mpconplane3d_->Evaluate(p,null,null,null,null,null);
    mpconline2d_->Evaluate(p,null,null,null,null,null);
    SynchronizeSumConstraint(p,actvalues_,"computed MPC value",numConstrID_,minConstrID_);
     
    constrainterr_->Update(1.0,*referencevalues_,-1.0,*actvalues_,0.0);
    return;
}


/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add scaled error of constraint to Lagrange multiplier.                 |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrMult(double factor)
{
  lagrMultVec_->Update(factor,*constrainterr_,1.0);
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add Lagrange increment to Lagrange multiplier.                         |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrMult(RCP<Epetra_Vector> vect)
{
  lagrMultVec_->Update(1.0,*vect,1.0);
  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
void ConstrManager::ComputeMonitorValues(RCP<Epetra_Vector> disp)
{
  vector<DRT::Condition*> monitcond(0);
  monitorvalues_->Scale(0.0);
  ParameterList p;
  actdisc_->SetState("displacement",disp);
 
  volmonitor3d_->Evaluate(p,null,null,null,null,null);
  SynchronizeSumConstraint(p, monitorvalues_,"computed volume",numMonitorID_,minMonitorID_);
  areamonitor3d_->Evaluate(p,null,null,null,null,null);
  areamonitor2d_->Evaluate(p,null,null,null,null,null);
  SynchronizeSumConstraint(p, monitorvalues_,"computed area",numMonitorID_,minMonitorID_);
  
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Print monitored values                                                 |
*-----------------------------------------------------------------------*/
void ConstrManager::PrintMonitorValues()
{
  for (int i = 0; i < numMonitorID_; ++i)
  {
    printf("Monitor value %2d: %10.5e (%5.2f%% of initial value)\n",i+minMonitorID_,(*monitorvalues_)[i],((*monitorvalues_)[i])*100/((*initialmonvalues_)[i]));
  }

  return;
}


/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraints by summing them up                                        |
 *----------------------------------------------------------------------*/
void ConstrManager::SynchronizeSumConstraint(ParameterList& params,
                                RCP<Epetra_Vector>& vect,
                                const char* resultstring, const int numID, const int minID)
{
  for (int i = 0; i < numID; ++i)
  {
    char valname[30];
    sprintf(valname,"%s %d",resultstring,i+minID);
    double currval=0.0;
    actdisc_->Comm().SumAll(&(params.get(valname,0.0)),&(currval),1);
    vect->SumIntoGlobalValues(1,&currval,&i);
  }
  return;
}


/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraint by finding minimum value                                   |
 *----------------------------------------------------------------------*/
void ConstrManager::SynchronizeMinConstraint(ParameterList& params,
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
