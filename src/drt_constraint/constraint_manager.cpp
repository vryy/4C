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
#include "../drt_lib/linalg_systemmatrix.H"
#include "iostream"
#include "../drt_fsi/fsi_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
ConstrManager::ConstrManager(DRT::Discretization& discr,
        RCP<Epetra_Vector> disp,
        ParameterList params):
actdisc_(rcp(&discr))
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
    actdisc_->EvaluateCondition(p,"VolumeConstraint_3D");
    havevolconstr_=true;
  }
  // Check for Area Constraints in 3D
  actdisc_->GetCondition("AreaConstraint_3D",constrcond);
  //Deal with area constraints in 3D
  if (constrcond.size())
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    actdisc_->EvaluateCondition(p,"AreaConstraint_3D");
    haveareaconstr3D_=true;
  }
  // Check for Area Constraints in 2D
  actdisc_->GetCondition("AreaConstraint_2D",constrcond);
  //Deal with area constraints
  if (constrcond.size())
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    actdisc_->EvaluateCondition(p,"AreaConstraint_2D");
    haveareaconstr2D_=true;
  }
  // Check for Nodal Constraints on planes in 3D
  //actdisc_->Print(cout);
  actdisc_->GetCondition("MPC_NodeOnPlane_3D",constrcond);
  //Deal with area constraints
  if (constrcond.size())
  {
    CreateDiscretizationFromCondition(constrcond,"ConstrDisc","CONSTRELE3");
    
//    
//    p.set("action","calc_nodeplanedist");
//    actdisc_->SetState("displacement",disp);
//    actdisc_->EvaluateCondition(p,"NodeOnPlaneConstraint_3D");
//    havenodeconstraint_=true;
  }
  //----------------------------------------------------
  //-----------include possible further constraints here
  //----------------------------------------------------
  haveconstraint_= haveareaconstr3D_||havevolconstr_||haveareaconstr2D_||havenodeconstraint_;
  if (haveconstraint_)
  {
    uzawaparam_=params.get<double>("uzawa parameter",1);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    ManageIDs(p,minConstrID_,maxConstrID_,numConstrID_);
    //initialize constraint Matrix
    constrMatrix_=rcp(new LINALG::SparseMatrix(*dofrowmap,numConstrID_,true,true));
    //build domainmap of constrMatrix
    constrmap_=rcp(new Epetra_Map(numConstrID_,0,actdisc_->Comm()));
    // sum up initial values
    initialvalues_=rcp(new Epetra_Vector(*constrmap_));
    initialvalues_->Scale(0.0);
    SynchronizeSumConstraint(p,initialvalues_,"computed volume",numConstrID_,minConstrID_);
    SynchronizeSumConstraint(p,initialvalues_,"computed area",numConstrID_,minConstrID_);
    SynchronizeSumConstraint(p,initialvalues_,"computed dist",numConstrID_,minConstrID_);
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
  havevolmonitor_=false;
  haveareamonitor3D_=false;
  haveareamonitor2D_=false;
  //Check for Volume Monitors
  actdisc_->GetCondition("VolumeMonitor_3D",constrcond);
  ParameterList p1;
  if (constrcond.size())
  {
    p1.set("action","calc_struct_constrvol");
    actdisc_->SetState("displacement",disp);
    actdisc_->EvaluateCondition(p1,"VolumeMonitor_3D");
    havevolmonitor_=true;
  }
  // Check for Area Monitor in 3D
  actdisc_->GetCondition("AreaMonitor_3D",constrcond);
  //Deal with area Monitors
  if (constrcond.size())
  {
    p1.set("action","calc_struct_monitarea");
    actdisc_->SetState("displacement",disp);
    actdisc_->EvaluateCondition(p1,"AreaMonitor_3D");
    haveareamonitor3D_=true;
  }
  // Check for Area Monitor in 2D
  actdisc_->GetCondition("AreaMonitor_2D",constrcond);
  //Deal with area Monitors
  if (constrcond.size())
  {
    p1.set("action","calc_struct_constrarea");
    actdisc_->SetState("displacement",disp);
    actdisc_->EvaluateCondition(p1,"AreaMonitor_2D");
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
  actvalues_->Scale(0.0);
  constrMatrix_->PutScalar(0.0);
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
  //Deal with volume constraints
  if (havevolconstr_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_volconstrstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",true);
    p.set("assemble vector 3",false);
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
    actdisc_->EvaluateCondition(p,stiff,constrMatrix_,fint,null,null,"VolumeConstraint_3D");
    SynchronizeSumConstraint(p,actvalues_,"computed volume",numConstrID_,minConstrID_);
  }
  //Deal with area constraints in 3D
  if (haveareaconstr3D_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_areaconstrstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",true);
    p.set("assemble vector 3",false);
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
    actdisc_->EvaluateCondition(p,stiff,constrMatrix_,fint,null,null,"AreaConstraint_3D");
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
  }
  //Deal with area constraints in 2D
  if (haveareaconstr2D_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_struct_areaconstrstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",true);
    p.set("assemble vector 3",false);
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
    actdisc_->EvaluateCondition(p,stiff,constrMatrix_,fint,null,null,"AreaConstraint_2D");
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
  }
  //Deal with area constraints in 2D
  if (havenodeconstraint_)
  {
    //Evaluate volume at predicted ENDpoint D_{n+1}
    // action for elements
    p.set("action","calc_areaconstrstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",true);
    p.set("assemble vector 3",false);
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
    actdisc_->EvaluateCondition(p,stiff,constrMatrix_,fint,null,null,"AreaConstraint_2D");
    SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
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
    actvalues_->Scale(0.0);
    ParameterList p;
    p.set("total time",time);
    actdisc_->SetState("displacement",disp);
    if(havevolconstr_)
    {
        p.set("action","calc_struct_constrvol");
        actdisc_->EvaluateCondition(p,"VolumeConstraint_3D");
        SynchronizeSumConstraint(p, actvalues_,"computed volume",numConstrID_,minConstrID_);
    }
    if(haveareaconstr3D_)
    {
        p.set("action","calc_struct_constrarea");
        actdisc_->EvaluateCondition(p,"AreaConstraint_3D");
        SynchronizeSumConstraint(p, actvalues_,"computed area",numConstrID_,minConstrID_);
    }
    if(haveareaconstr2D_)
    {
      p.set("action","calc_struct_constrarea");
      actdisc_->EvaluateCondition(p,"AreaConstraint_2D");
      SynchronizeSumConstraint(p, actvalues_,"computed area",numConstrID_,minConstrID_);
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
  monitorvalues_->Scale(0.0);
  ParameterList p;
  actdisc_->SetState("displacement",disp);
  if(havevolmonitor_)
  {
    p.set("action","calc_struct_constrvol");
    actdisc_->EvaluateCondition(p,"VolumeMonitor_3D");
    SynchronizeSumConstraint(p, monitorvalues_,"computed volume",numMonitorID_,minMonitorID_);
  }
  if(haveareamonitor3D_)
  {
    p.set("action","calc_struct_monitarea");
    actdisc_->EvaluateCondition(p,"AreaMonitor_3D");
    SynchronizeSumConstraint(p, monitorvalues_,"computed area",numMonitorID_,minMonitorID_);
  }
  if(haveareamonitor2D_)
  {
    p.set("action","calc_struct_constrarea");
    actdisc_->EvaluateCondition(p,"AreaMonitor_2D");
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
 |(private)                                                 tk 11/07    |
 |small subroutine for initialization purposes to determine ID's        |
 *----------------------------------------------------------------------*/
void ConstrManager::ManageIDs(ParameterList& params, int& minID, int& maxID, int& numID)
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


RCP<DRT::Discretization> ConstrManager::CreateDiscretizationFromCondition
        (   vector< DRT::Condition* >      constrcondvec,
            const string&             discret_name,
            const string&             element_name)
{
  RCP<Epetra_Comm> com = rcp(actdisc_->Comm().Clone());
  
  constraintdis_ = rcp(new DRT::Discretization(discret_name,com));

  if (!actdisc_->Filled()) actdisc_->FillComplete(true,true,true,false);

  const int myrank = constraintdis_->Comm().MyPID();
  
 if(constrcondvec.size()==0)
      dserror("number of multi point constraint conditions = 0 --> cannot create constraint discretization");
  
  // vector with boundary ele id's
  vector<int> egid;
  //egid.reserve(actdisc_->NumMyRowElements());
//
  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* actnoderowmap = actdisc_->NodeRowMap();
//  
  // Loop all conditions is constrcondvec
  for (unsigned int j=0;j<constrcondvec.size();j++)
  {

    vector<int> ngid=*(constrcondvec[j]->Nodes());
    const int numnodes=ngid.size();
    // We sort the global node ids according to the definition of the boundary condition
    SortConstraintNodes(ngid, constrcondvec[j]);
     
    remove_copy_if(&ngid[0], &ngid[0]+numnodes,
                     inserter(rownodeset, rownodeset.begin()),
                     not1(FSI::UTILS::MyGID(actnoderowmap)));
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
        constraintdis_->AddNode(rcp(new DRT::Node(gid, standardnode->X(), myrank)));
      }
    }
    //build unique node row map
    vector<int> boundarynoderowvec(rownodeset.begin(), rownodeset.end());
    rownodeset.clear();
    RCP<Epetra_Map> constraintnoderowmap = rcp(new Epetra_Map(-1,
                                                               boundarynoderowvec.size(),
                                                               &boundarynoderowvec[0],
                                                               0,
                                                               constraintdis_->Comm()));
    boundarynoderowvec.clear();

    //build overlapping node column map
    vector<int> constraintnodecolvec(colnodeset.begin(), colnodeset.end());
    colnodeset.clear();
    RCP<Epetra_Map> constraintnodecolmap = rcp(new Epetra_Map(-1,
                                                               constraintnodecolvec.size(),
                                                               &constraintnodecolvec[0],
                                                               0,
                                                               constraintdis_->Comm()));
    constraintnodecolvec.clear();
    
    
    int myrank = actdisc_->Comm().MyPID();
      // create an element with the same global element id
    if (actdisc_->gNode( ngid[0] )->Owner() == myrank )
    {
      RCP<DRT::Element> constraintele = DRT::UTILS::Factory(element_name, j, myrank);
  
      // set the same global node ids to the ale element
      constraintele->SetNodeIds(ngid.size(), &(ngid[0]));
  
      // add constraint element
      
      constraintdis_->AddElement(constraintele);

    }
    // now care about the parallel distribution and ghosting.
    // So far every processor only knows about his nodes
        
    constraintdis_->ExportColumnNodes(*constraintnodecolmap);
  
    RefCountPtr< Epetra_Map > constraintelerowmap;
    RefCountPtr< Epetra_Map > constraintelecolmap;
  
    // now we have all elements in a linear map roweles
    // build resonable maps for elements from the
    // already valid and final node maps
    // note that nothing is actually redistributed in here
    constraintdis_->BuildElementRowColumn(*constraintnoderowmap, *constraintnodecolmap, constraintelerowmap, constraintelecolmap);
  
    // we can now export elements to resonable row element distribution
    constraintdis_->ExportRowElements(*constraintelerowmap);
  
    // export to the column map / create ghosting of elements
    constraintdis_->ExportColumnElements(*constraintelecolmap);
  
  }
  // Now we are done. :)
 
  constraintdis_->FillComplete(true,true,true,false);
  constraintdis_->Print(cout);  
  return null;
}

void ConstrManager::SortConstraintNodes(
    vector<int>& nodeids,
    DRT::Condition*      cond)
{

  // get this condition's nodes
  vector<int> temp=nodeids;
  const vector<int>*    NodeID1  = cond->Get<vector<int> >("NodeID 1");
  const vector<int>*    NodeID2  = cond->Get<vector<int> >("NodeID 2");
  const vector<int>*    NodeID3  = cond->Get<vector<int> >("NodeID 3");
  const vector<int>*    NodeID4  = cond->Get<vector<int> >("NodeID 4");
  nodeids[0]=temp[(*NodeID1)[0]-1];
  nodeids[1]=temp[(*NodeID2)[0]-1];
  nodeids[2]=temp[(*NodeID3)[0]-1];
  nodeids[3]=temp[(*NodeID4)[0]-1];
  return;
}

#endif
