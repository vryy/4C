/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2

 *----------------------------------------------------------------------*/


#include "4C_constraint_multipointconstraint3penalty.hpp"

#include "4C_constraint_element3.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_dofset_transparent.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_utils_function_of_time.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint3Penalty::MPConstraint3Penalty(
    Teuchos::RCP<DRT::Discretization> discr,  ///< discretization constraint lives on
    const std::string& CondName               ///< Name of condition to create constraint from
    )
    : MPConstraint(discr, CondName)
{
  if (constrcond_.size())
  {
    // control the constraint by absolute or relative values
    for (auto* conditer : constrcond_)
    {
      const int condID = conditer->parameters().Get<int>("ConditionID");
      penalties_[condID] = conditer->parameters().Get<double>("penalty");
      const std::string type = conditer->parameters().Get<std::string>("control");
      if (type == "abs")
        absconstraint_[condID] = true;
      else
      {
        absconstraint_[condID] = false;
      }
    }

    int startID = 0;
    constraintdis_ = create_discretization_from_condition(
        actdisc_, constrcond_, "ConstrDisc", "CONSTRELE3", startID);

    std::map<int, Teuchos::RCP<DRT::Discretization>>::iterator discriter;
    for (discriter = constraintdis_.begin(); discriter != constraintdis_.end(); discriter++)
    {
      Teuchos::RCP<Epetra_Map> newcolnodemap =
          CORE::REBALANCE::ComputeNodeColMap(actdisc_, discriter->second);
      actdisc_->Redistribute(*(actdisc_->NodeRowMap()), *newcolnodemap);
      Teuchos::RCP<CORE::Dofsets::DofSet> newdofset =
          Teuchos::rcp(new CORE::Dofsets::TransparentDofSet(actdisc_));
      (discriter->second)->ReplaceDofSet(newdofset);
      newdofset = Teuchos::null;
      (discriter->second)->fill_complete();
    }

    int nummyele = 0;
    int numele = eletocond_id_.size();
    if (!actdisc_->Comm().MyPID())
    {
      nummyele = numele;
    }
    // initialize maps and importer
    errormap_ = Teuchos::rcp(new Epetra_Map(numele, nummyele, 0, actdisc_->Comm()));
    rederrormap_ = CORE::LINALG::AllreduceEMap(*errormap_);
    errorexport_ = Teuchos::rcp(new Epetra_Export(*rederrormap_, *errormap_));
    errorimport_ = Teuchos::rcp(new Epetra_Import(*rederrormap_, *errormap_));
    acterror_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    initerror_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
  }
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::Initialize(const double& time)
{
  for (auto* cond : constrcond_)
  {
    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond->parameters().Get<int>("ConditionID");

    // if current time (at) is larger than activation time of the condition, activate it
    if ((inittimes_.find(condID)->second < time) && (!activecons_.find(condID)->second))
    {
      activecons_.find(condID)->second = true;
      if (actdisc_->Comm().MyPID() == 0)
      {
        std::cout << "Encountered another active condition (Id = " << condID
                  << ")  for restart time t = " << time << std::endl;
      }
    }
  }
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::Initialize(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  FOUR_C_THROW("method not used for penalty formulation!");
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::Initialize(Teuchos::ParameterList& params)
{
  const double time = params.get("total time", -1.0);

  for (auto* cond : constrcond_)
  {
    int condID = cond->parameters().Get<int>("ConditionID");
    // control absolute values
    switch (Type())
    {
      case mpcnodeonplane3d:
      case mpcnormalcomp3d:
        params.set("action", "calc_MPC_state");
        break;
      case none:
        return;
      default:
        FOUR_C_THROW("Constraint/monitor is not an multi point constraint!");
    }

    evaluate_error(constraintdis_.find(condID)->second, params, initerror_, true);

    activecons_.find(condID)->second = true;
    if (actdisc_->Comm().MyPID() == 0)
    {
      std::cout << "Encountered a new active condition (Id = " << condID
                << ")  at time t = " << time << std::endl;
    }
    //    std::cout << "initial error "<< *initerror_<<std::endl;
  }
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::Evaluate(Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  switch (Type())
  {
    case mpcnodeonplane3d:
    case mpcnormalcomp3d:
      params.set("action", "calc_MPC_state");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Constraint/monitor is not an multi point constraint!");
  }

  acterror_->PutScalar(0.0);
  std::map<int, Teuchos::RCP<DRT::Discretization>>::iterator discriter;
  for (discriter = constraintdis_.begin(); discriter != constraintdis_.end(); discriter++)

    evaluate_error(discriter->second, params, acterror_);

  //    std::cout << "current error "<< *acterror_<<std::endl;

  switch (Type())
  {
    case mpcnodeonplane3d:
    case mpcnormalcomp3d:
      params.set("action", "calc_MPC_stiff");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Constraint/monitor is not an multi point constraint!");
  }
  for (discriter = constraintdis_.begin(); discriter != constraintdis_.end(); discriter++)
    evaluate_constraint(discriter->second, params, systemmatrix1, systemmatrix2, systemvector1,
        systemvector2, systemvector3);

  return;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::map<int, Teuchos::RCP<DRT::Discretization>>
CONSTRAINTS::MPConstraint3Penalty::create_discretization_from_condition(
    Teuchos::RCP<DRT::Discretization> actdisc,
    std::vector<CORE::Conditions::Condition*> constrcondvec, const std::string& discret_name,
    const std::string& element_name, int& startID)
{
  // start with empty map
  std::map<int, Teuchos::RCP<DRT::Discretization>> newdiscmap;

  if (!actdisc->Filled())
  {
    actdisc->fill_complete();
  }

  if (constrcondvec.size() == 0)
    FOUR_C_THROW(
        "number of multi point constraint conditions = 0 --> cannot create constraint "
        "discretization");

  // Loop all conditions in constrcondvec and build discretization for any condition ID

  int index = 0;  // counter for the index of condition in vector
  std::vector<CORE::Conditions::Condition*>::iterator conditer;
  for (conditer = constrcondvec.begin(); conditer != constrcondvec.end(); conditer++)
  {
    // initialize a new discretization
    Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(actdisc->Comm().Clone());
    Teuchos::RCP<DRT::Discretization> newdis =
        Teuchos::rcp(new DRT::Discretization(discret_name, com, actdisc->n_dim()));
    const int myrank = newdis->Comm().MyPID();
    std::set<int> rownodeset;
    std::set<int> colnodeset;
    const Epetra_Map* actnoderowmap = actdisc->NodeRowMap();
    // get node IDs, this vector will only contain FREE nodes in the end
    std::vector<int> ngid = *((*conditer)->GetNodes());
    std::vector<int> defnv;
    switch (Type())
    {
      case mpcnodeonplane3d:
      {
        // take three nodes defining plane as specified by user and put them into a set
        const auto& defnvp = (*conditer)->parameters().Get<std::vector<int>>("planeNodes");
        defnv = defnvp;
      }
      break;
      case mpcnormalcomp3d:
      {
        // take master node
        const int defn = (*conditer)->parameters().Get<int>("masterNode");
        defnv.push_back(defn);
      }
      break;
      default:
        FOUR_C_THROW("not good!");
    }
    std::set<int> defns(defnv.begin(), defnv.end());
    std::set<int>::iterator nsit;
    // safe gids of definition nodes in a vector
    std::vector<int> defnodeIDs;

    int counter = 1;  // counter is used to keep track of deleted node ids from the vector, input
                      // starts with 1

    for (nsit = defns.begin(); nsit != defns.end(); ++nsit)
    {
      defnodeIDs.push_back(ngid.at((*nsit) - counter));
      ngid.erase(ngid.begin() + (*nsit) - counter);
      counter++;
    }

    unsigned int nodeiter;
    // loop over all free nodes of condition
    for (nodeiter = 0; nodeiter < ngid.size(); nodeiter++)
    {
      std::vector<int> ngid_ele = defnodeIDs;
      ngid_ele.push_back(ngid[nodeiter]);
      const int numnodes = ngid_ele.size();
      remove_copy_if(ngid_ele.data(), ngid_ele.data() + numnodes,
          inserter(rownodeset, rownodeset.begin()),
          std::not_fn(CORE::Conditions::MyGID(actnoderowmap)));
      // copy node ids specified in condition to colnodeset
      copy(ngid_ele.data(), ngid_ele.data() + numnodes, inserter(colnodeset, colnodeset.begin()));

      // construct constraint nodes, which use the same global id as the standard nodes
      for (int i = 0; i < actnoderowmap->NumMyElements(); ++i)
      {
        const int gid = actnoderowmap->GID(i);
        if (rownodeset.find(gid) != rownodeset.end())
        {
          const CORE::Nodes::Node* standardnode = actdisc->lRowNode(i);
          newdis->AddNode(Teuchos::rcp(new CORE::Nodes::Node(gid, standardnode->X(), myrank)));
        }
      }

      if (myrank == 0)
      {
        Teuchos::RCP<CORE::Elements::Element> constraintele =
            CORE::COMM::Factory(element_name, "Polynomial", nodeiter + startID, myrank);
        // set the same global node ids to the ale element
        constraintele->SetNodeIds(ngid_ele.size(), ngid_ele.data());
        // add constraint element
        newdis->add_element(constraintele);
      }
      // save the connection between element and condition
      eletocond_id_[nodeiter + startID] = (*conditer)->parameters().Get<int>("ConditionID");
      eletocondvecindex_[nodeiter + startID] = index;
    }
    // adjust starting ID for next condition, in this case nodeiter=ngid.size(), hence the counter
    // is larger than the ID
    // of the last element
    startID += nodeiter;

    // now care about the parallel distribution and ghosting.
    // So far every processor only knows about his nodes

    // build unique node row map
    std::vector<int> boundarynoderowvec(rownodeset.begin(), rownodeset.end());
    rownodeset.clear();
    Teuchos::RCP<Epetra_Map> constraintnoderowmap = Teuchos::rcp(new Epetra_Map(
        -1, boundarynoderowvec.size(), boundarynoderowvec.data(), 0, newdis->Comm()));
    boundarynoderowvec.clear();

    // build overlapping node column map
    std::vector<int> constraintnodecolvec(colnodeset.begin(), colnodeset.end());
    colnodeset.clear();
    Teuchos::RCP<Epetra_Map> constraintnodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, constraintnodecolvec.size(), constraintnodecolvec.data(), 0, newdis->Comm()));

    constraintnodecolvec.clear();
    newdis->Redistribute(*constraintnoderowmap, *constraintnodecolmap);
    // put new discretization into the map
    newdiscmap[(*conditer)->parameters().Get<int>("ConditionID")] = newdis;
    // increase counter
    index++;
  }

  startID--;  // set counter back to ID of the last element
  return newdiscmap;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::evaluate_constraint(Teuchos::RCP<DRT::Discretization> disc,
    Teuchos::ParameterList& params, Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  if (!(disc->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!(disc->HaveDofs())) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  // define element matrices and vectors
  CORE::LINALG::SerialDenseMatrix elematrix1;
  CORE::LINALG::SerialDenseMatrix elematrix2;
  CORE::LINALG::SerialDenseVector elevector1;
  CORE::LINALG::SerialDenseVector elevector2;
  CORE::LINALG::SerialDenseVector elevector3;

  const double time = params.get("total time", -1.0);
  const int numcolele = disc->NumMyColElements();

  // get values from time integrator to scale matrices with
  double scStiff = params.get("scaleStiffEntries", 1.0);

  // loop over column elements
  for (int i = 0; i < numcolele; ++i)
  {
    // some useful data for computation
    CORE::Elements::Element* actele = disc->lColElement(i);
    int eid = actele->Id();
    int condID = eletocond_id_.find(eid)->second;
    CORE::Conditions::Condition* cond = constrcond_[eletocondvecindex_.find(eid)->second];
    params.set<Teuchos::RCP<CORE::Conditions::Condition>>("condition", Teuchos::rcp(cond, false));

    // computation only if time is larger or equal than initialization time for constraint
    if (inittimes_.find(condID)->second <= time)
    {
      // initialize if it is the first time condition is evaluated
      if (activecons_.find(condID)->second == false)
      {
        const std::string action = params.get<std::string>("action");
        Teuchos::RCP<Epetra_Vector> displast = params.get<Teuchos::RCP<Epetra_Vector>>("old disp");
        SetConstrState("displacement", displast);
        // last converged step is used reference
        Initialize(params);
        Teuchos::RCP<Epetra_Vector> disp = params.get<Teuchos::RCP<Epetra_Vector>>("new disp");
        SetConstrState("displacement", disp);
        params.set("action", action);
      }

      // get element location vector, dirichlet flags and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      actele->LocationVector(*disc, lm, lmowner, lmstride);
      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();
      elematrix1.shape(eledim, eledim);
      elevector1.size(eledim);
      elevector3.size(1);
      params.set("ConditionID", eid);

      // call the element evaluate method
      int err = actele->Evaluate(
          params, *disc, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("Proc %d: Element %d returned err=%d", disc->Comm().MyPID(), eid, err);

      // loadcurve business
      const auto* curve = cond->parameters().GetIf<int>("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve);
      double curvefac = 1.0;
      bool usetime = true;
      if (time < 0.0) usetime = false;
      if (curvenum >= 0 && usetime)
        curvefac = GLOBAL::Problem::Instance()
                       ->FunctionById<CORE::UTILS::FunctionOfTime>(curvenum)
                       .Evaluate(time);


      double diff = (curvefac * (*initerror_)[eid] - (*acterror_)[eid]);
      elematrix1.scale(diff);
      for (int i = 0; i < eledim; i++)
        for (int j = 0; j < eledim; j++) elematrix1(i, j) += elevector1(i) * elevector1(j);
      elematrix1.scale(2 * scStiff * penalties_[condID]);

      systemmatrix1->Assemble(eid, lmstride, elematrix1, lm, lmowner);
      elevector1.scale(2. * penalties_[condID] * diff);
      CORE::LINALG::Assemble(*systemvector1, elevector1, lm, lmowner);
    }
  }
}  // end of evaluate_condition

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint3Penalty::evaluate_error(Teuchos::RCP<DRT::Discretization> disc,
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector, bool init)
{
  if (!(disc->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!(disc->HaveDofs())) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  // define element matrices and vectors
  CORE::LINALG::SerialDenseMatrix elematrix1;
  CORE::LINALG::SerialDenseMatrix elematrix2;
  CORE::LINALG::SerialDenseVector elevector1;
  CORE::LINALG::SerialDenseVector elevector2;
  CORE::LINALG::SerialDenseVector elevector3;

  // loop over column elements
  const double time = params.get("total time", -1.0);
  const int numcolele = disc->NumMyColElements();
  for (int i = 0; i < numcolele; ++i)
  {
    // some useful data for computation
    CORE::Elements::Element* actele = disc->lColElement(i);
    int eid = actele->Id();
    int condID = eletocond_id_.find(eid)->second;
    CORE::Conditions::Condition* cond = constrcond_[eletocondvecindex_.find(eid)->second];
    params.set<Teuchos::RCP<CORE::Conditions::Condition>>("condition", Teuchos::rcp(cond, false));

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->LocationVector(*disc, lm, lmowner, lmstride);
    elevector3.size(1);
    params.set("ConditionID", eid);

    if (absconstraint_.find(condID)->second && init)
    {
      elevector3[0] = cond->parameters().Get<double>("amplitude");
    }
    else
    {
      // call the element evaluate method
      int err = actele->Evaluate(
          params, *disc, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("Proc %d: Element %d returned err=%d", disc->Comm().MyPID(), eid, err);
    }

    // assembly
    std::vector<int> constrlm;
    std::vector<int> constrowner;
    constrlm.push_back(eid);
    constrowner.push_back(actele->Owner());
    CORE::LINALG::Assemble(*systemvector, elevector3, constrlm, constrowner);

    activecons_.find(condID)->second = true;

    if (actdisc_->Comm().MyPID() == 0 && (!(activecons_.find(condID)->second)))
    {
      std::cout << "Encountered a new active penalty mp condition (Id = " << condID
                << ")  at time t = " << time << std::endl;
    }
  }

  Teuchos::RCP<Epetra_Vector> acterrdist = Teuchos::rcp(new Epetra_Vector(*errormap_));
  acterrdist->Export(*systemvector, *errorexport_, Add);
  systemvector->Import(*acterrdist, *errorimport_, Insert);
  return;
}  // end of evaluate_error

FOUR_C_NAMESPACE_CLOSE
