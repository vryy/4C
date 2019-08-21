/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with constraints living on boundaries
\level 2

\maintainer Matthias Mayr

*----------------------------------------------------------------------*/


#include "monitor.H"
#include <iostream>
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::Monitor::Monitor(Teuchos::RCP<DRT::Discretization> discr, const std::string& conditionname,
    int& minID, int& maxID)
    : actdisc_(discr)
{
  actdisc_->GetCondition(conditionname, moncond_);
  if (moncond_.size())
  {
    montype_ = GetMoniType(conditionname);
    for (unsigned int i = 0; i < moncond_.size(); i++)
    {
      // moncond_[i]->Print(std::cout);
      int condID = (*(moncond_[i]->Get<std::vector<int>>("ConditionID")))[0];
      if (condID > maxID)
      {
        maxID = condID;
      }
      if (condID < minID)
      {
        minID = condID;
      }
    }
  }
  else
  {
    montype_ = none;
  }
}


/*-----------------------------------------------------------------------*
|(private)                                                       tk 07/08|
*-----------------------------------------------------------------------*/
UTILS::Monitor::MoniType UTILS::Monitor::GetMoniType(const std::string& name)
{
  if (name == "VolumeMonitor_3D")
    return volmonitor3d;
  else if (name == "AreaMonitor_3D")
    return areamonitor3d;
  else if (name == "AreaMonitor_2D")
    return areamonitor2d;
  return none;
}


/*-----------------------------------------------------------------------*
|(public)                                                        tk 07/08|
|Evaluate Monitors, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void UTILS::Monitor::Evaluate(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  switch (montype_)
  {
    case volmonitor3d:
      params.set("action", "calc_struct_constrvol");
      break;
    case areamonitor3d:
      params.set("action", "calc_struct_monitarea");
      break;
    case areamonitor2d:
      params.set("action", "calc_struct_constrarea");
      break;
    case none:
      return;
    default:
      dserror("Unknown monitor type to be evaluated in Monitor class!");
  }
  EvaluateMonitor(params, systemvector);
  return;
}


/*-----------------------------------------------------------------------*
 |(private)                                                     tk 08/08 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void UTILS::Monitor::EvaluateMonitor(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < moncond_.size(); ++i)
  {
    DRT::Condition& cond = *(moncond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    const std::vector<int>* CondIDVec = cond.Get<std::vector<int>>("ConditionID");
    const int condID = (*CondIDVec)[0];
    const int offsetID = params.get("OffsetID", 0);
    params.set<Teuchos::RCP<DRT::Condition>>("condition", Teuchos::rcp(&cond, false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_, lm, lmowner, lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      elevector3.Size(1);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(
          params, *actdisc_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly
      std::vector<int> constrlm;
      std::vector<int> constrowner;
      constrlm.push_back(condID - offsetID);
      constrowner.push_back(curr->second->Owner());
      LINALG::Assemble(*systemvector, elevector3, constrlm, constrowner);
    }
  }
  return;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Monitor::SetState(const std::string& state,  ///< name of state to set
    Teuchos::RCP<Epetra_Vector> V                        ///< values to set
)
{
  actdisc_->SetState(state, V);
}
