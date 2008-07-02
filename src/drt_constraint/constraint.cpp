/*!----------------------------------------------------------------------
\file constraint.cpp

\brief Basic constraint class, dealing with constraints living on boundaries
<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraint.H"
#include "constraint_element.H"
#include "mpcdofset.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
Constraint::Constraint(RCP<DRT::Discretization> discr,
        const string& conditionname,
        int& minID,
        int& maxID):
actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,constrcond_);
  if (constrcond_.size())
  {
    constrtype_=GetConstrType(conditionname);
    for (unsigned int i=0; i<constrcond_.size();i++)
    {
      int condID=(*(constrcond_[i]->Get<vector<int> >("ConditionID")))[0];
      if (condID>maxID)
        maxID=condID;
      if (condID<minID)
        minID=condID;
    }
  }
  else
  {
    constrtype_=none;
  }
  
}


/*-----------------------------------------------------------------------*
|(private)                                                       tk 07/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
Constraint::ConstrType Constraint::GetConstrType(const string& name)
{
  if (name=="VolumeConstraint_3D")
    return volconstr3d;
  else if (name=="AreaConstraint_3D")
    return areaconstr3d;
  else if (name=="AreaConstraint_2D")
    return areaconstr2d;
  else if (name=="VolumeMonitor_3D")
    return volmonitor3d;
  else if (name=="AreaMonitor_3D")
    return areamonitor3d;
  else if (name=="AreaMonitor_2D")
    return areamonitor2d;
  else if (name=="MPC_NodeOnPlane_3D")
    return mpcnodeonplane3d;
  else if (name=="MPC_NodeOnLine_2D")
    return mpcnodeonline2d;
  
  return none;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 07/08|
|Evaluate Constraints, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void Constraint::Evaluate(
    ParameterList&        params,
    RCP<LINALG::SparseOperator> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3)
{
  
  if ((systemmatrix1==Teuchos::null)&&(systemmatrix2==Teuchos::null)&&
      (systemvector1==Teuchos::null)&&(systemvector2==Teuchos::null))
  {
    switch (constrtype_)
    {
      case volconstr3d: 
        params.set("action","calc_struct_constrvol");
      break;
      case areaconstr3d:
        params.set("action","calc_struct_constrarea");
      break;
      case areaconstr2d:
        params.set("action","calc_struct_constrarea");
      break;
      case volmonitor3d:
        params.set("action","calc_struct_constrvol");
      break;
      case areamonitor3d:
        params.set("action","calc_struct_monitarea");
      break;
      case areamonitor2d:
        params.set("action","calc_struct_constrarea");
      break;
      case none:
        return;
      default:
        dserror("Unknown constraint/monitor type to be evaluated in Constraint class!");
    }
  }
  else switch (constrtype_)
  {
    case volconstr3d: 
      params.set("action","calc_struct_volconstrstiff");
    break;
    case areaconstr3d:
      params.set("action","calc_struct_areaconstrstiff");
    break;
    case areaconstr2d:
      params.set("action","calc_struct_areaconstrstiff");
    break;
    case none:
      return;  
    default:
      dserror("Wrong constraint/monitor type to evaluate systemvector!");
  }
  EvaluateConstraint(actdisc_,params,systemmatrix1,systemmatrix2,systemvector1,systemvector2,systemvector3);
  return;
}

/*-----------------------------------------------------------------------*
 |(private)                                                     tk 07/08 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void Constraint::EvaluateConstraint(RCP<DRT::Discretization> disc,
    ParameterList&        params,
    RCP<LINALG::SparseOperator> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3)
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
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
      DRT::Condition& cond = *(constrcond_[i]);

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
        elevector3.Size(systemvector3->MyLength());

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*disc,lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly
        int eid = curr->second->Id();
        if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,lm,lmowner);
        if (assemblemat2)
        {
          //assemble to rectangular matrix. The colum corresponds to the constraint ID
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
            constrowner.push_back(curr->second->Owner());
          }
          LINALG::Assemble(*systemvector3,elevector3,constrlm,constrowner);
        }
      }
  }
  return;
} // end of EvaluateCondition

#endif
