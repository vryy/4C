/*!----------------------------------------------------------------------
\file multipointconstraint.cpp

\brief Basic constraint class, dealing with multi point constraints
<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET


#include "multipointconstraint.H"
#include "mpcdofset.H"
#include "constraint_element.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::MPConstraint::MPConstraint(RCP<DRT::Discretization> discr,
        const string& conditionname,
        int& minID,
        int& maxID)
: UTILS::Constraint
  (
    discr,
    conditionname,
    minID,
    maxID
  )
{
  if (constrcond_.size())
  {
    constraintdis_=CreateDiscretizationFromCondition(actdisc_,constrcond_,"ConstrDisc","CONSTRELE");
    ReplaceNumDof(actdisc_,constraintdis_);
    RCP<Epetra_Map> newcolnodemap = ComputeNodeColMap(actdisc_, constraintdis_);
    actdisc_->Redistribute(*(actdisc_->NodeRowMap()), *newcolnodemap);
    RCP<DRT::DofSet> newdofset = rcp(new MPCDofSet(actdisc_));
    constraintdis_->ReplaceDofSet(newdofset);
    newdofset = null;
    constraintdis_->FillComplete();
  }
  return;
}


/*-----------------------------------------------------------------------*
|(public)                                                        tk 07/08|
|Evaluate Constraints, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void UTILS::MPConstraint::Evaluate(
    ParameterList&        params,
    RCP<LINALG::SparseOperator> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3,
    bool                  init)
{ 
  // in case init is set to true we want to set systemvector1 to the amplitudes defined 
  // in the input file
  if (init)
  {
    //allocate vectors for amplitudes and IDs
    vector<double> amplit(constrcond_.size());
    vector<int> IDs(constrcond_.size());
    
    //read data of the input files
    for (unsigned int i=0;i<constrcond_.size();i++)
    {
      const vector<double>*    MPCampl  = constrcond_[i]->Get<vector<double> >("Amplitude");
      const vector<int>*    MPCcondID  = constrcond_[i]->Get<vector<int> >("ConditionID");
      amplit[i]=(*MPCampl)[0];
      const int mid=params.get("MinID",0);
      IDs[i]=(*MPCcondID)[0]-mid;
    }
    // replace systemvector1 by the given amplitude values
    // systemvector1 is supposed to be the vector with initial values of the constraints
    systemvector1->ReplaceGlobalValues(amplit.size(),&(amplit[0]),&(IDs[0]));
  }
  else 
  {
    switch (Type())
    {
      case mpcnodeonplane3d: 
      case mpcnodeonline2d:
        params.set("action","calc_MPC_stiff");
      break;
      case none:
        return;
      default:
        dserror("Constraint/monitor is not an multi point constraint!");
    }
    EvaluateConstraint(constraintdis_,params,systemmatrix1,systemmatrix2,systemvector1,systemvector2,systemvector3);
  }
  return;
}

/*------------------------------------------------------------------------*
 |(private)                                                   tk 04/08    |
 |subroutine creating a new discretization containing constraint elements |
 *------------------------------------------------------------------------*/
RCP<DRT::Discretization> UTILS::MPConstraint::CreateDiscretizationFromCondition
(  
  RCP<DRT::Discretization> actdisc,
  vector< DRT::Condition* >      constrcondvec,
  const string&             discret_name,
  const string&             element_name
)
{
  RCP<Epetra_Comm> com = rcp(actdisc->Comm().Clone());

  RCP<DRT::Discretization> newdis = rcp(new DRT::Discretization(discret_name,com));

  if (!actdisc->Filled())
  {
    actdisc->FillComplete();
  }

  const int myrank = newdis->Comm().MyPID();

  if(constrcondvec.size()==0)
      dserror("number of multi point constraint conditions = 0 --> cannot create constraint discretization");

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* actnoderowmap = actdisc->NodeRowMap();

  // Loop all conditions in constrcondvec
  for (unsigned int j=0;j<constrcondvec.size();j++)
  {
    vector<int> ngid=*(constrcondvec[j]->Nodes());
    const int numnodes=ngid.size();
    // We sort the global node ids according to the definition of the boundary condition
    ReorderConstraintNodes(ngid, constrcondvec[j]);

    remove_copy_if(&ngid[0], &ngid[0]+numnodes,
                     inserter(rownodeset, rownodeset.begin()),
                     not1(DRT::UTILS::MyGID(actnoderowmap)));
    // copy node ids specified in condition to colnodeset
    copy(&ngid[0], &ngid[0]+numnodes,
          inserter(colnodeset, colnodeset.begin()));

    // construct boundary nodes, which use the same global id as the cutter nodes
    for (int i=0; i<actnoderowmap->NumMyElements(); ++i)
    {
      const int gid = actnoderowmap->GID(i);
      if (rownodeset.find(gid)!=rownodeset.end())
      {
        const DRT::Node* standardnode = actdisc->lRowNode(i);
        newdis->AddNode(rcp(new DRT::Node(gid, standardnode->X(), myrank)));
      }
    }

    if (myrank == 0)
    {
      RCP<DRT::Element> constraintele = DRT::UTILS::Factory(element_name,"Polynomial", j, myrank);
      // set the same global node ids to the ale element
      constraintele->SetNodeIds(ngid.size(), &(ngid[0]));
      
      // add constraint element
      newdis->AddElement(constraintele);

    }
    // now care about the parallel distribution and ghosting.
    // So far every processor only knows about his nodes
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

  DRT::UTILS::RedistributeWithNewNodalDistribution(*newdis,*constraintnoderowmap,*constraintnodecolmap);
  
  return newdis;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 04/08    |
 |reorder MPC nodes based on condition input                            |
 *----------------------------------------------------------------------*/
void UTILS::MPConstraint::ReorderConstraintNodes
(
  vector<int>& nodeids,
  const DRT::Condition* cond
)
{
  // get this condition's nodes
  vector<int> temp=nodeids;
  if (nodeids.size()==4)
  {
    const vector<int>*    constrNode  = cond->Get<vector<int> >("ConstrNode");
    nodeids[(*constrNode)[0]-1]=temp[3];
    nodeids[3]=temp[(*constrNode)[0]-1];
  }
  else if (nodeids.size()==3)
  {
    const vector<int>*    constrNode1  = cond->Get<vector<int> >("ConstrNode 1");
    const vector<int>*    constrNode2  = cond->Get<vector<int> >("ConstrNode 2");
    const vector<int>*    constrNode3  = cond->Get<vector<int> >("ConstrNode 3");
    nodeids[0]=temp[(*constrNode1)[0]-1];
    nodeids[1]=temp[(*constrNode2)[0]-1];
    nodeids[2]=temp[(*constrNode3)[0]-1];
  }
  else
  {
    dserror("strange number of nodes for an MPC! Should be 3 (2D) or 4 (3D)");
  }
  return;
}


/*----------------------------------------------------------------------*
 |(private)                                                   tk 05/08  |
 |replace numdofs in elements of constraint discretization              |
 *----------------------------------------------------------------------*/
void UTILS::MPConstraint::ReplaceNumDof
(
  const RCP<DRT::Discretization> sourcedis,
  const RCP<DRT::Discretization> constraintdis
) const
{
  // find typical numdof of basis discretization (may not work for XFEM)
  const DRT::Element* actele = sourcedis->lColElement(0);
  const DRT::Node*const* nodes = actele->Nodes();
  const int mpc_numdof = sourcedis->NumDof(nodes[0]);

  // change numdof for all constraint elements
  const int numcolele = constraintdis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::ELEMENTS::ConstraintElement* mpcele = dynamic_cast<DRT::ELEMENTS::ConstraintElement*>(constraintdis->lColElement(i));
    mpcele->SetNumDofPerNode(mpc_numdof);
  }
  constraintdis->FillComplete();
  return;
 }


/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |recompute nodecolmap of standard discretization to include constrained|
 |nodes as ghosted nodes                                                |
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> UTILS::MPConstraint::ComputeNodeColMap(
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

/*-----------------------------------------------------------------------*
 |(private)                                                     tk 07/08 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void UTILS::MPConstraint::EvaluateConstraint(RCP<DRT::Discretization> disc,
    ParameterList&        params,
    RCP<LINALG::SparseOperator> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3)
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
    if (assemblemat1) elematrix1.Shape(eledim,eledim);
    if (assemblemat2) elematrix2.Shape(eledim,eledim);
    if (assemblevec1) elevector1.Size(eledim);
    if (assemblevec2) elevector2.Size(eledim);        
    if (assemblevec3) elevector3.Size(systemvector3->MyLength());
    elevector3.Size(systemvector3->MyLength());
    DRT::Condition& cond = *(constrcond_[actele->Id()]);
    const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
    int condID=(*CondIDVec)[0];
    params.set("ConditionID",condID);
    params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));
    // call the element evaluate method
    int err = actele->Evaluate(params,*disc,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",disc->Comm().MyPID(),actele->Id(),err);

    int eid = actele->Id();
    if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,lm,lmowner);
    if (assemblemat2)
    {
      int minID=params.get("MinID",0);
      vector<int> colvec(1);
      colvec[0]=condID-minID;
      systemmatrix2->Assemble(eid,elevector2,lm,lmowner,colvec);
    }
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) 
    {
      vector<int> constrlm;
      vector<int> constrowner;
      for (int i=0; i<elevector3.Length();i++)
      {
        constrlm.push_back(i);
        constrowner.push_back(actele->Owner());
      }
      LINALG::Assemble(*systemvector3,elevector3,constrlm,constrowner);
    }
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
} // end of EvaluateCondition

#endif
