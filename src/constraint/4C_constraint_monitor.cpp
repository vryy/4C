/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with constraints living on boundaries
\level 2


*----------------------------------------------------------------------*/


#include "4C_constraint_monitor.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::Monitor::Monitor(Teuchos::RCP<Core::FE::Discretization> discr,
    const std::string& conditionname, int& minID, int& maxID)
    : actdisc_(discr)
{
  actdisc_->GetCondition(conditionname, moncond_);
  if (moncond_.size())
  {
    montype_ = get_moni_type(conditionname);
    for (auto& i : moncond_)
    {
      int condID = i->parameters().get<int>("ConditionID");

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
CONSTRAINTS::Monitor::MoniType CONSTRAINTS::Monitor::get_moni_type(const std::string& name)
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
void CONSTRAINTS::Monitor::evaluate(
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
      FOUR_C_THROW("Unknown monitor type to be evaluated in Monitor class!");
  }
  evaluate_monitor(params, systemvector);
}


/*-----------------------------------------------------------------------*
 |(private)                                                     tk 08/08 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void CONSTRAINTS::Monitor::evaluate_monitor(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  if (!(actdisc_->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!actdisc_->HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (auto* cond : moncond_)
  {
    // Get ConditionID of current condition if defined and write value in parameterlist
    const int condID = cond->parameters().get<int>("ConditionID");
    const int offsetID = params.get("OffsetID", 0);
    params.set<Teuchos::RCP<Core::Conditions::Condition>>("condition", Teuchos::rcp(cond, false));

    // define element matrices and vectors
    Core::LinAlg::SerialDenseMatrix elematrix1;
    Core::LinAlg::SerialDenseMatrix elematrix2;
    Core::LinAlg::SerialDenseVector elevector1;
    Core::LinAlg::SerialDenseVector elevector2;
    Core::LinAlg::SerialDenseVector elevector3;

    std::map<int, Teuchos::RCP<Core::Elements::Element>>& geom = cond->Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator curr;
    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_, lm, lmowner, lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      elevector3.size(1);

      // call the element specific evaluate method
      int err = curr->second->evaluate(
          params, *actdisc_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("error while evaluating elements");

      // assembly
      std::vector<int> constrlm;
      std::vector<int> constrowner;
      constrlm.push_back(condID - offsetID);
      constrowner.push_back(curr->second->Owner());
      Core::LinAlg::Assemble(*systemvector, elevector3, constrlm, constrowner);
    }
  }
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::Monitor::set_state(const std::string& state,  ///< name of state to set
    Teuchos::RCP<Epetra_Vector> V                               ///< values to set
)
{
  actdisc_->set_state(state, V);
}

FOUR_C_NAMESPACE_CLOSE
