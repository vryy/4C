/*!----------------------------------------------------------------------
\file drt_constraint_manager.cpp

\brief Class controlling constraint and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraint_manager.H"
#include "constraint_element3.H"
#include "mpcdofset.H"
#include "iostream"
#include "../drt_adapter/adapter_utils.H"
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
  haveareaconstr3D_=false;
  haveareaconstr2D_=false;
  havevolconstr_=false;
  havenodeconstraint_=false;
  //Check for volume constraints
  vector<DRT::Condition*> constrcond(0);
  actdisc_->GetCondition("VolumeConstraint_3D",constrcond);
  // Keep ParameterList p alive during initialization, so global information
  // over IDs as well as element results stored here can be used after all
  // constraints are evaluated
  ParameterList p;
  //Deal with volume constraints
  if (constrcond.size())
  {
    p.set("action","calc_struct_constrvol");
    actdisc_->SetState("displacement",disp);
//    actdisc_->EvaluateCondition(p,"VolumeConstraint_3D");
    EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
    havevolconstr_=true;
  }
  // Check for Area Constraints in 3D
  actdisc_->GetCondition("AreaConstraint_3D",constrcond);
  //Deal with area constraints in 3D
  if (constrcond.size())
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
    haveareaconstr3D_=true;
  }
  // Check for Area Constraints in 2D
  actdisc_->GetCondition("AreaConstraint_2D",constrcond);
  //Deal with area constraints
  if (constrcond.size())
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
    haveareaconstr2D_=true;
  }
  // Check for Multi Point Constraint Node on planes in 3D
  actdisc_->GetCondition("MPC_NodeOnPlane_3D",constrcond);
  //Initialize Vectors to contain IDs and amplitudes
  vector<double> MPCamplitudes(constrcond.size());
  vector<int> MPCcondIDs(constrcond.size());
  // Deal with MPC
  if (constrcond.size())
  {
    // Construct special constraint discretization consisting of constraint elements
    constraintdis_=CreateDiscretizationFromCondition(constrcond,"ConstrDisc","CONSTRELE3");
    // find typical numdof of basis discretization (may not work for XFEM)
    const DRT::Element* actele = actdisc_->lColElement(0);
    const DRT::Node*const* nodes = actele->Nodes();
    const int mpc_numdof = actdisc_->NumDof(nodes[0]);
    
    // change numdof for all constraint elements
    const int numcolele = constraintdis_->NumMyColElements();
    for (int i=0; i<numcolele; ++i)
    {
      DRT::ELEMENTS::ConstraintElement* mpcele = dynamic_cast<DRT::ELEMENTS::ConstraintElement*>(constraintdis_->lColElement(i));
      mpcele->SetNumDofPerNode(mpc_numdof);
    }
    constraintdis_->FillComplete();
    // now reconstruct the extended colmap
    RCP<Epetra_Map> newcolnodemap = ComputeNodeColMap(actdisc_, constraintdis_);
    //Redistribute underlying discretization to have ghosting in the same way as in the
    //constraint discretization
    actdisc_->Redistribute(*(actdisc_->NodeRowMap()), *newcolnodemap);
    // Change dof's in constraint discretization to fit underlying discretization
    RCP<DRT::DofSet> newdofset = rcp(new MPCDofSet(actdisc_));
    constraintdis_->ReplaceDofSet(newdofset);
    newdofset = null;
    constraintdis_->FillComplete();
    //Fill MPC Vectors
    SetupMPC(MPCamplitudes,MPCcondIDs,constrcond,p);
    havenodeconstraint_=true;
  }
  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  haveconstraint_= haveareaconstr3D_||havevolconstr_||haveareaconstr2D_||havenodeconstraint_;
  if (haveconstraint_)
  {
    uzawaparam_=params.get<double>("uzawa parameter",1);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    ManageIDs(p,minConstrID_,maxConstrID_,numConstrID_,MPCcondIDs);
    //initialize constraint Matrix
    constrMatrix_=rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,true,true));
    //build domainmap of constrMatrix
    constrmap_=rcp(new Epetra_Map(numConstrID_,0,actdisc_->Comm()));
    // sum up initial values
    initialvalues_=rcp(new Epetra_Vector(*constrmap_));
    initialvalues_->Scale(0.0);
    //Set initial values to computed volumes and areas and to amplitudes of MPC
    SynchronizeSumConstraint(p,initialvalues_,"computed volume",numConstrID_,minConstrID_);
    SynchronizeSumConstraint(p,initialvalues_,"computed area",numConstrID_,minConstrID_);

    initialvalues_->ReplaceGlobalValues(MPCamplitudes.size(),&(MPCamplitudes[0]),&(MPCcondIDs[0]));
    //Initialize Lagrange Multiplicators, reference values and errors
    actdisc_->ClearState();
    referencevalues_=rcp(new Epetra_Vector(*constrmap_));
    actvalues_=rcp(new Epetra_Vector(*constrmap_));
    actvalues_->Scale(0.0);
    constrainterr_=rcp(new Epetra_Vector(*constrmap_));
    lagrMultVec_=rcp(new Epetra_Vector(*constrmap_));
    lagrMultInc_=rcp(new Epetra_Vector(*constrmap_));
    lagrMultVec_->Scale(0.0);
    lagrMultInc_->Scale(0.0);
    fact_=rcp(new Epetra_Vector(*constrmap_));
  }
  //-----------------------------Monitor Conditions!
  havevolmonitor_=false;
  haveareamonitor3D_=false;
  haveareamonitor2D_=false;
  vector<DRT::Condition*> monitcond(0);
  //Check for Volume Monitors
  actdisc_->GetCondition("VolumeMonitor_3D",monitcond);
  ParameterList p1;
  if (monitcond.size())
  {
    p1.set("action","calc_struct_constrvol");
    actdisc_->SetState("displacement",disp);
    EvaluateCondition(actdisc_,p1,null,null,null,null,null,monitcond);
    havevolmonitor_=true;
  }
  // Check for Area Monitor in 3D
  actdisc_->GetCondition("AreaMonitor_3D",monitcond);
  //Deal with area Monitors
  if (monitcond.size())
  {
    p1.set("action","calc_struct_monitarea");
    actdisc_->SetState("displacement",disp);
    EvaluateCondition(actdisc_,p1,null,null,null,null,null,monitcond);
    haveareamonitor3D_=true;
  }
  // Check for Area Monitor in 2D
  actdisc_->GetCondition("AreaMonitor_2D",monitcond);
  //Deal with area Monitors
  if (monitcond.size())
  {
    p1.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    EvaluateCondition(actdisc_,p1,null,null,null,null,null,monitcond);
    haveareamonitor2D_=true;
  }
  //----------------------------------------------------
  //--------------include possible further monitors here
  //----------------------------------------------------
  havemonitor_= haveareamonitor3D_||havevolmonitor_||haveareamonitor2D_;
  if (havemonitor_)
  {
    ManageIDs(p1,minMonitorID_,maxMonitorID_,numMonitorID_);
    // sum up initial values
    monitormap_=rcp(new Epetra_Map(numMonitorID_,0,actdisc_->Comm()));
    monitorvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_=rcp(new Epetra_Vector(*monitormap_));
    initialmonvalues_->Scale(0.0);
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
  //Deal with volume constraints
  if (havevolconstr_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_volconstrstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("MinID",minConstrID_);
    p.set("MaxID",maxConstrID_);
    p.set("NumberofID",numConstrID_);
    // Convert Epetra_Vector constaining lagrange multipliers to an completely
    // redundant Epetra_vector since every element with the constraint condition needs them
    RCP<Epetra_Map> reducedmap = LINALG::AllreduceEMap(*constrmap_);
    RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*reducedmap));
    LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
    p.set("LagrMultVector",lagrMultVecDense);
    actdisc_->ClearState();
    actdisc_->SetState("displacement",disp);
    actdisc_->GetCondition("VolumeConstraint_3D",constrcond);
    EvaluateCondition(actdisc_,p,stiff,constrMatrix_,fint,null,null,constrcond);
    SynchronizeSumConstraint(p,actvalues_,"computed volume",numConstrID_,minConstrID_);
  }
  //Deal with area constraints in 3D
  if (haveareaconstr3D_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_areaconstrstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("MinID",minConstrID_);
    p.set("MaxID",maxConstrID_);
    p.set("NumberofID",numConstrID_);
    // Convert Epetra_Vector constaining lagrange multipliers to an completely
    // redundant Epetra_vector since every element with the constraint condition needs them
    RCP<Epetra_Map> reducedmap = LINALG::AllreduceEMap(*constrmap_);
    RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*reducedmap));
    LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
    p.set("LagrMultVector",lagrMultVecDense);
    actdisc_->ClearState();
    actdisc_->SetState("displacement",disp);
    actdisc_->GetCondition("AreaConstraint_3D",constrcond);
    EvaluateCondition(actdisc_,p,stiff,constrMatrix_,fint,null,null,constrcond);
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
  }
  //Deal with area constraints in 2D
  if (haveareaconstr2D_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_areaconstrstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("MinID",minConstrID_);
    p.set("MaxID",maxConstrID_);
    p.set("NumberofID",numConstrID_);
    // Convert Epetra_Vector constaining lagrange multipliers to an completely
    // redundant Epetra_vector since every element with the constraint condition needs them
    RCP<Epetra_Map> reducedmap = LINALG::AllreduceEMap(*constrmap_);
    RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*reducedmap));
    LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
    p.set("LagrMultVector",lagrMultVecDense);
    actdisc_->ClearState();
    actdisc_->SetState("displacement",disp);
    actdisc_->GetCondition("AreaConstraint_2D",constrcond);
    EvaluateCondition(actdisc_,p,stiff,constrMatrix_,fint,null,null,constrcond);
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
  }
  //Deal with MPC
  if (havenodeconstraint_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_MPC_stiff");
    // other parameters that are needed by the constraint elements
    p.set("total time",time);
    p.set("MinID",minConstrID_);
    p.set("MaxID",maxConstrID_);
    p.set("NumberofID",numConstrID_);
    // Convert Epetra_Vector constaining lagrange multipliers to an completely
    // redundant Epetra_vector since every element with the constraint condition needs them
    RCP<Epetra_Map> reducedmap = LINALG::AllreduceEMap(*constrmap_);
    RCP<Epetra_Vector> lagrMultVecDense = rcp(new Epetra_Vector(*reducedmap));
    LINALG::Export(*lagrMultVec_,*lagrMultVecDense);
    p.set("LagrMultVector",lagrMultVecDense);
    constraintdis_->ClearState();
    constraintdis_->SetState("displacement",disp);
    //conditions to evaluate
    vector<DRT::Condition*> constrcond(0);
    actdisc_->GetCondition("MPC_NodeOnPlane_3D",constrcond);
    Evaluate(constraintdis_,p,stiff,constrMatrix_,fint,null,null,constrcond);
    SynchronizeSumConstraint(p,actvalues_,"computed normal distance",numConstrID_,minConstrID_);
  }

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
    if(havevolconstr_)
    {
        p.set("action","calc_struct_constrvol");
        actdisc_->GetCondition("VolumeConstraint_3D",constrcond);
        EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
        SynchronizeSumConstraint(p, actvalues_,"computed volume",numConstrID_,minConstrID_);
    }
    if(haveareaconstr3D_)
    {
        p.set("action","calc_struct_constrarea");
        actdisc_->GetCondition("AreaConstraint_3D",constrcond);
        EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
        SynchronizeSumConstraint(p, actvalues_,"computed area",numConstrID_,minConstrID_);
    }
    if(haveareaconstr2D_)
    {
      p.set("action","calc_struct_constrarea");
      actdisc_->GetCondition("AreaConstraint_2D",constrcond);
      EvaluateCondition(actdisc_,p,null,null,null,null,null,constrcond);
      SynchronizeSumConstraint(p, actvalues_,"computed area",numConstrID_,minConstrID_);
    }
    if(havenodeconstraint_)
    {
      p.set("action","calc_MPC_stiff");
      actdisc_->GetCondition("MPC_NodeOnPlane_3D",constrcond);
      Evaluate(constraintdis_,p,null,null,null,null,null,constrcond);
      SynchronizeSumConstraint(p,actvalues_,"computed normal distance",numConstrID_,minConstrID_);
    }
    constrainterr_->Update(1.0,*referencevalues_,-1.0,*actvalues_,0.0);
    return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Update Lagrange increment according to linear uzawa algorithm.          |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrIncr(double factor, Epetra_Vector vect)
{
  lagrMultInc_->Update(factor,vect,factor,*constrainterr_,1.0);
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
void ConstrManager::UpdateLagrMult()
{
  lagrMultVec_->Update(1.0,*lagrMultInc_,1.0);
  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute matrix-vector-product of matrix of constraint conditions with   |
|lagrange increment as needed for uzawa algorithm                        |
*------------------------------------------------------------------------*/
void ConstrManager::ComputeConstrTimesLagrIncr(RCP<Epetra_Vector> dotprod)
{
  dotprod->PutScalar(0.0);
  constrMatrix_->Multiply(false,*lagrMultInc_,*dotprod) ;
  dotprod->Scale(-1.0);
  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute matrix-vector-product of matrix of constraint conditions with   |
|displacement increment as needed to compute residual of uzawa algorithm |
*-----------------------------------------------------------------------*/
void ConstrManager::ComputeConstrTimesDisi(Epetra_Vector disi, RCP<Epetra_Vector> dotprod)
{

  dotprod->PutScalar(0.0);
  constrMatrix_->Multiply(true,disi,*dotprod) ;
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
  if(havevolmonitor_)
  {
    p.set("action","calc_struct_constrvol");
    actdisc_->GetCondition("VolumeConstraint_3D",monitcond);
    EvaluateCondition(actdisc_,p,null,null,null,null,null,monitcond);
    SynchronizeSumConstraint(p, monitorvalues_,"computed volume",numMonitorID_,minMonitorID_);
  }
  if(haveareamonitor3D_)
  {
    p.set("action","calc_struct_monitarea");
    actdisc_->GetCondition("AreaMonitor_3D",monitcond);
    EvaluateCondition(actdisc_,p,null,null,null,null,null,monitcond);
    SynchronizeSumConstraint(p, monitorvalues_,"computed area",numMonitorID_,minMonitorID_);
  }
  if(haveareamonitor2D_)
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->GetCondition("AreaMonitor_2D",monitcond);
    EvaluateCondition(actdisc_,p,null,null,null,null,null,monitcond);
    SynchronizeSumConstraint(p, monitorvalues_,"computed area",numMonitorID_,minMonitorID_);
  }
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
 |(private)                                                 tk 04/08    |
 |small subroutine for initialization purposes to determine ID's        |
 *----------------------------------------------------------------------*/
void ConstrManager::ManageIDs(ParameterList& params,
    int& minID,
    int& maxID,
    int& numID,
    vector<int>& MPCIds)
{
  actdisc_->Comm().MaxAll(&(params.get("MaxID",0)),&maxID,1);
  actdisc_->Comm().MinAll(&(params.get("MinID",100000)),&minID,1);
  numID=1+maxID-minID;
  for (unsigned int var = 0; var < MPCIds.size(); ++var)
  {
    MPCIds[var]-=minID;
  }
  return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine for initialization purposes to determine ID's        |
 *----------------------------------------------------------------------*/
void ConstrManager::ManageIDs(ParameterList& params,
    int& minID,
    int& maxID,
    int& numID)
{
  actdisc_->Comm().MaxAll(&(params.get("MaxID",0)),&maxID,1);
  actdisc_->Comm().MinAll(&(params.get("MinID",100000)),&minID,1);
  numID=1+maxID-minID;

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

/*------------------------------------------------------------------------*
 |(private)                                                   tk 04/08    |
 |subroutine creating a new discretization containing constraint elements |
 *------------------------------------------------------------------------*/
RCP<DRT::Discretization> ConstrManager::CreateDiscretizationFromCondition
        (   vector< DRT::Condition* >      constrcondvec,
            const string&             discret_name,
            const string&             element_name)
{
  RCP<Epetra_Comm> com = rcp(actdisc_->Comm().Clone());
  
  RCP<DRT::Discretization> newdis = rcp(new DRT::Discretization(discret_name,com));

  if (!actdisc_->Filled())
  {
    actdisc_->FillComplete();
  }

  const int myrank = newdis->Comm().MyPID();

  if(constrcondvec.size()==0)
      dserror("number of multi point constraint conditions = 0 --> cannot create constraint discretization");

  // vector with boundary ele id's
  //vector<int> egid;

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* actnoderowmap = actdisc_->NodeRowMap();
//
  // Loop all conditions in constrcondvec
  for (unsigned int j=0;j<constrcondvec.size();j++)
  {

    vector<int> ngid=*(constrcondvec[j]->Nodes());
    const int numnodes=ngid.size();
    // We sort the global node ids according to the definition of the boundary condition
    ReorderConstraintNodes(ngid, constrcondvec[j]);

    remove_copy_if(&ngid[0], &ngid[0]+numnodes,
                     inserter(rownodeset, rownodeset.begin()),
                     not1(ADAPTER::UTILS::MyGID(actnoderowmap)));
    // copy node ids specified in condition to colnodeset
    copy(&ngid[0], &ngid[0]+numnodes,
          inserter(colnodeset, colnodeset.begin()));

    // construct boundary nodes, which use the same global id as the cutter nodes
    for (int i=0; i<actnoderowmap->NumMyElements(); ++i)
    {
      const int gid = actnoderowmap->GID(i);
      if (rownodeset.find(gid)!=rownodeset.end())
      {
        const DRT::Node* standardnode = actdisc_->lRowNode(i);
        newdis->AddNode(rcp(new DRT::Node(gid, standardnode->X(), myrank)));
      }
    }
    //build unique node row map
    vector<int> boundarynoderowvec(rownodeset.begin(), rownodeset.end());
    rownodeset.clear();
    RCP<Epetra_Map> constraintnoderowmap = rcp(new Epetra_Map(-1,
                                                               boundarynoderowvec.size(),
                                                               &boundarynoderowvec[0],
                                                               0,
                                                               newdis->Comm()));
    boundarynoderowvec.clear();

    //build overlapping node column map
    vector<int> constraintnodecolvec(colnodeset.begin(), colnodeset.end());
    colnodeset.clear();
    RCP<Epetra_Map> constraintnodecolmap = rcp(new Epetra_Map(-1,
                                                               constraintnodecolvec.size(),
                                                               &constraintnodecolvec[0],
                                                               0,
                                                               newdis->Comm()));
    constraintnodecolvec.clear();


    int myrank = actdisc_->Comm().MyPID();

    if (myrank == 0)
    {
      RCP<DRT::Element> constraintele = DRT::UTILS::Factory(element_name, j, myrank);
      // set the same global node ids to the ale element
      constraintele->SetNodeIds(ngid.size(), &(ngid[0]));

      // add constraint element
      newdis->AddElement(constraintele);

    }
    // now care about the parallel distribution and ghosting.
    // So far every processor only knows about his nodes

    newdis->ExportColumnNodes(*constraintnodecolmap);

    RefCountPtr< Epetra_Map > constraintelerowmap;
    RefCountPtr< Epetra_Map > constraintelecolmap;

    // now we have all elements in a linear map roweles
    // build resonable maps for elements from the
    // already valid and final node maps
    // note that nothing is actually redistributed in here
    newdis->BuildElementRowColumn(*constraintnoderowmap, *constraintnodecolmap, constraintelerowmap, constraintelecolmap);

    // we can now export elements to resonable row element distribution
    newdis->ExportRowElements(*constraintelerowmap);

    // export to the column map / create ghosting of elements
    newdis->ExportColumnElements(*constraintelecolmap);

  }
  newdis->FillComplete();
  return newdis;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 04/08    |
 |reorder MPC nodes based on condition input                            |
 *----------------------------------------------------------------------*/
void ConstrManager::ReorderConstraintNodes(
    vector<int>& nodeids,
    const DRT::Condition*      cond)
{
  // get this condition's nodes
  vector<int> temp=nodeids;
  const vector<int>*    constrNode  = cond->Get<vector<int> >("ConstrNode");
  nodeids[(*constrNode)[0]-1]=temp[3];
  nodeids[3]=temp[(*constrNode)[0]-1];
  return;
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |recompute nodecolmap of standard discretization to include constrained|
 |nodes as ghosted nodes                                                |
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> ConstrManager::ComputeNodeColMap(
        const RCP<DRT::Discretization> sourcedis,
        const RCP<DRT::Discretization> constraintdis
        ) const
{
    const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

    vector<int> mycolnodes(oldcolnodemap->NumMyElements());
    oldcolnodemap->MyGlobalElements (&mycolnodes[0]);
    for (int inode = 0; inode != constraintdis->NumMyColNodes(); ++inode)
    {
        const DRT::Node* newnode = constraintdis->lColNode(inode);
        const int gid = newnode->Id();
        if (!(sourcedis->HaveGlobalNode(gid)))
        {
            mycolnodes.push_back(gid);
        }
    }

    // now reconstruct the extended colmap
    RCP<Epetra_Map> newcolnodemap = rcp(new Epetra_Map(-1,
                                       mycolnodes.size(),
                                       &mycolnodes[0],
                                       0,
                                       sourcedis->Comm()));
    return newcolnodemap;
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |fill vectors with MPC Ids and amplitudes and update parameter list    |
 *----------------------------------------------------------------------*/
void ConstrManager::SetupMPC(vector<double>& amplit,
    vector<int>& IDs,
    const vector<DRT::Condition*>& constrcond,
    Teuchos::ParameterList& p)
{
  for (unsigned int i=0;i<constrcond.size();i++)
  {
    const vector<double>*    MPCampl  = constrcond[i]->Get<vector<double> >("Amplitude");
    const vector<int>*    MPCcondID  = constrcond[i]->Get<vector<int> >("ConditionID");
    amplit[i]=(*MPCampl)[0];
    IDs[i]=(*MPCcondID)[0];
    // element write condition IDs (min and max) to parameter list. Here we have to do it ourselves
    // since we do not call any element action
    const int maxID=p.get("MaxID",0);
    const int minID=p.get("MinID",1000000);
    if (maxID<(*MPCcondID)[0])
    {
      p.set("MaxID",(*MPCcondID)[0]);
    }
    if (minID>(*MPCcondID)[0])
    {
      p.set("MinID",(*MPCcondID)[0]);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |Evaluate method, calling element evaluates of a discretization and    |
 |assembing results based on given conditions                           |
 *----------------------------------------------------------------------*/
void ConstrManager::Evaluate( RCP<DRT::Discretization> disc,
                              ParameterList&        params,
                              RCP<LINALG::SparseOperator> systemmatrix1,
                              RCP<LINALG::SparseOperator> systemmatrix2,
                              RCP<Epetra_Vector>    systemvector1,
                              RCP<Epetra_Vector>    systemvector2,
                              RCP<Epetra_Vector>    systemvector3,
                              vector<DRT::Condition*>& constrcond)
{
  if (!(disc->Filled())) dserror("FillComplete() was not called");
  if (!(disc->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblemat2 = systemmatrix2!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over column elements
  const int numcolele = disc->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = disc->lColElement(i);

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(*disc,lm,lmowner);

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim);
    elematrix2.Shape(eledim,eledim);
    elevector1.Size(eledim);
    elevector2.Size(eledim);
    elevector3.Size(eledim);
    DRT::Condition& cond = *(constrcond[actele->Id()]);
    const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
    int condID=(*CondIDVec)[0];
    params.set("ConditionID",condID);
    params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));
    // call the element evaluate method
    int err = actele->Evaluate(params,*disc,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",disc->Comm().MyPID(),actele->Id(),err);

    if (assemblemat1) systemmatrix1->Assemble(elematrix1,lm,lmowner);
    if (assemblemat2)
    {
      int minID=params.get("MinID",0);
      vector<int> colvec(1);
      colvec[0]=condID-minID;
      systemmatrix2->Assemble(elevector2,lm,lmowner,colvec);
    }
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);
    
    const vector<int>*    curve  = cond.Get<vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    bool usetime = true;
    const double time = params.get("total time",-1.0);
    if (time<0.0) usetime = false;
    if (curvenum>=0 && usetime)
      curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

    // Get ConditionID of current condition if defined and write value in parameterlist
    char factorname[30];
    sprintf(factorname,"LoadCurveFactor %d",condID);
    params.set(factorname,curvefac);
  }

  return;
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |Evaluate method, calling element evaluates of a condition and         |
 |assembing results based on this conditions                            |
 *----------------------------------------------------------------------*/
void ConstrManager::EvaluateCondition(RCP<DRT::Discretization> disc,
    ParameterList&        params,
    RCP<LINALG::SparseOperator> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3,
    vector<DRT::Condition*>& constrcond)
{
  if (!(disc->Filled())) dserror("FillComplete() was not called");
  if (!disc->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond.size(); ++i) 
  {
      DRT::Condition& cond = *(constrcond[i]);

      // Evaluate Loadcurve if defined. Put current load factor in parameterlist
      const vector<int>*    curve  = cond.Get<vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

      // Get ConditionID of current condition if defined and write value in parameterlist
      
      const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
      int condID=(*CondIDVec)[0];
      params.set("ConditionID",condID);
      char factorname[30];
      sprintf(factorname,"LoadCurveFactor %d",condID);
      params.set(factorname,curvefac);
      
      params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // if (geom.empty()) dserror("evaluation of condition with empty geometry");
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*disc,lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.Shape(eledim,eledim);
        elematrix2.Shape(eledim,eledim);
        elevector1.Size(eledim);
        elevector2.Size(eledim);
        elevector3.Size(eledim);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*disc,lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly
        if (assemblemat1) systemmatrix1->Assemble(elematrix1,lm,lmowner);
        if (assemblemat2)
        {
          int minID=params.get("MinID",0);
          vector<int> colvec(1);
          colvec[0]=condID-minID;
          systemmatrix2->Assemble(elevector2,lm,lmowner,colvec);
        }
        if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
        if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
        if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);
      }    
  } 
  return;
} // end of EvaluateCondition

#endif
