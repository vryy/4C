/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/



#include "4C_constraint_multipointconstraint2.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_transparent.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_utils_function_of_time.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint2::MPConstraint2(Teuchos::RCP<Core::FE::Discretization> discr,
    const std::string& conditionname, int& minID, int& maxID)
    : MPConstraint(discr, conditionname, minID, maxID)
{
  if (constrcond_.size())
  {
    int dummy = 0;
    // create constraint discretization and store it with label 0, within the map
    constraintdis_ = create_discretization_from_condition(
        actdisc_, constrcond_, "ConstrDisc", "CONSTRELE2", dummy);
    Teuchos::RCP<Epetra_Map> newcolnodemap =
        Core::Rebalance::ComputeNodeColMap(actdisc_, constraintdis_.find(0)->second);
    actdisc_->redistribute(*(actdisc_->node_row_map()), *newcolnodemap);
    Teuchos::RCP<Core::DOFSets::DofSet> newdofset =
        Teuchos::rcp(new Core::DOFSets::TransparentDofSet(actdisc_));
    (constraintdis_.find(0)->second)->replace_dof_set(newdofset);
    newdofset = Teuchos::null;
    (constraintdis_.find(0)->second)->fill_complete();
  }
}

/*------------------------------------------------------------------------*
|(public)                                                       tk 08/08  |
|Initialization routine activates conditions (restart)                    |
*------------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint2::initialize(const double& time)
{
  for (auto* cond : constrcond_)
  {
    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond->parameters().get<int>("ConditionID");

    // if current time (at) is larger than activation time of the condition, activate it
    if ((inittimes_.find(condID)->second < time) && (!activecons_.find(condID)->second))
    {
      activecons_.find(condID)->second = true;
      if (actdisc_->get_comm().MyPID() == 0)
      {
        std::cout << "Encountered another active condition (Id = " << condID
                  << ")  for restart time t = " << time << std::endl;
      }
    }
  }
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 07/08|
|Evaluate Constraints, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint2::initialize(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  const double time = params.get("total time", -1.0);
  // in case init is set to true we want to set systemvector1 to the amplitudes defined
  // in the input file
  // allocate vectors for amplitudes and IDs


  std::vector<double> amplit(constrcond_.size());
  std::vector<int> IDs(constrcond_.size());
  // read data of the input files
  for (unsigned int i = 0; i < constrcond_.size(); i++)
  {
    Core::Conditions::Condition& cond = *(constrcond_[i]);
    int condID = cond.parameters().get<int>("ConditionID");
    if (inittimes_.find(condID)->second <= time)
    {
      const int MPCcondID = constrcond_[i]->parameters().get<int>("ConditionID");
      amplit[i] = constrcond_[i]->parameters().get<double>("amplitude");
      const int mid = params.get("OffsetID", 0);
      IDs[i] = MPCcondID - mid;
      // remember next time, that this condition is already initialized, i.e. active
      activecons_.find(condID)->second = true;
      if (actdisc_->get_comm().MyPID() == 0)
      {
        std::cout << "Encountered a new active condition (Id = " << condID
                  << ")  at time t = " << time << std::endl;
      }
    }
  }
  // replace systemvector by the given amplitude values
  // systemvector is supposed to be the vector with initial values of the constraints
  if (actdisc_->get_comm().MyPID() == 0)
  {
    systemvector->ReplaceGlobalValues(amplit.size(), amplit.data(), IDs.data());
  }
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 07/08|
|Evaluate Constraints, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint2::evaluate(Teuchos::ParameterList& params,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  switch (type())
  {
    case mpcnodeonline2d:
      params.set("action", "calc_MPC_stiff");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Constraint/monitor is not an multi point constraint!");
  }
  evaluate_constraint(constraintdis_.find(0)->second, params, systemmatrix1, systemmatrix2,
      systemvector1, systemvector2, systemvector3);
}

/*------------------------------------------------------------------------*
 |(private)                                                   tk 04/08    |
 |subroutine creating a new discretization containing constraint elements |
 *------------------------------------------------------------------------*/
std::map<int, Teuchos::RCP<Core::FE::Discretization>>
CONSTRAINTS::MPConstraint2::create_discretization_from_condition(
    Teuchos::RCP<Core::FE::Discretization> actdisc,
    std::vector<Core::Conditions::Condition*> constrcondvec, const std::string& discret_name,
    const std::string& element_name, int& startID)
{
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(actdisc->get_comm().Clone());

  Teuchos::RCP<Core::FE::Discretization> newdis =
      Teuchos::rcp(new Core::FE::Discretization(discret_name, com, actdisc->n_dim()));

  if (!actdisc->filled())
  {
    actdisc->fill_complete();
  }

  const int myrank = newdis->get_comm().MyPID();

  if (constrcondvec.size() == 0)
    FOUR_C_THROW(
        "number of multi point constraint conditions = 0 --> cannot create constraint "
        "discretization");

  std::set<int> rownodeset;
  std::set<int> colnodeset;
  const Epetra_Map* actnoderowmap = actdisc->node_row_map();

  // Loop all conditions in constrcondvec
  for (unsigned int j = 0; j < constrcondvec.size(); j++)
  {
    std::vector<int> ngid = *(constrcondvec[j]->get_nodes());
    const int numnodes = ngid.size();
    // We sort the global node ids according to the definition of the boundary condition
    reorder_constraint_nodes(ngid, constrcondvec[j]);

    remove_copy_if(ngid.data(), ngid.data() + numnodes, inserter(rownodeset, rownodeset.begin()),
        std::not_fn(Core::Conditions::MyGID(actnoderowmap)));
    // copy node ids specified in condition to colnodeset
    copy(ngid.data(), ngid.data() + numnodes, inserter(colnodeset, colnodeset.begin()));

    // construct boundary nodes, which use the same global id as the cutter nodes
    for (int i = 0; i < actnoderowmap->NumMyElements(); ++i)
    {
      const int gid = actnoderowmap->GID(i);
      if (rownodeset.find(gid) != rownodeset.end())
      {
        const Core::Nodes::Node* standardnode = actdisc->l_row_node(i);
        newdis->add_node(Teuchos::rcp(new Core::Nodes::Node(gid, standardnode->x(), myrank)));
      }
    }

    if (myrank == 0)
    {
      Teuchos::RCP<Core::Elements::Element> constraintele =
          Core::Communication::Factory(element_name, "Polynomial", j, myrank);
      // set the same global node ids to the ale element
      constraintele->set_node_ids(ngid.size(), ngid.data());

      // add constraint element
      newdis->add_element(constraintele);
    }
    // now care about the parallel distribution and ghosting.
    // So far every processor only knows about his nodes
  }

  // build unique node row map
  std::vector<int> boundarynoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  Teuchos::RCP<Epetra_Map> constraintnoderowmap = Teuchos::rcp(new Epetra_Map(
      -1, boundarynoderowvec.size(), boundarynoderowvec.data(), 0, newdis->get_comm()));
  boundarynoderowvec.clear();

  // build overlapping node column map
  std::vector<int> constraintnodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  Teuchos::RCP<Epetra_Map> constraintnodecolmap = Teuchos::rcp(new Epetra_Map(
      -1, constraintnodecolvec.size(), constraintnodecolvec.data(), 0, newdis->get_comm()));

  constraintnodecolvec.clear();

  newdis->redistribute(*constraintnoderowmap, *constraintnodecolmap);

  std::map<int, Teuchos::RCP<Core::FE::Discretization>> newdismap;
  newdismap[startID] = newdis;
  return newdismap;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 04/08    |
 |reorder MPC nodes based on condition input                            |
 *----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint2::reorder_constraint_nodes(
    std::vector<int>& nodeids, const Core::Conditions::Condition* cond)
{
  // get this condition's nodes
  std::vector<int> temp = nodeids;
  if (nodeids.size() == 3)
  {
    nodeids[0] = temp[cond->parameters().get<int>("constrNode 1") - 1];
    nodeids[1] = temp[cond->parameters().get<int>("constrNode 2") - 1];
    nodeids[2] = temp[cond->parameters().get<int>("constrNode 3") - 1];
  }
  else
  {
    FOUR_C_THROW("strange number of nodes for an MPC! Should be 3 in 2D.");
  }
}

/*-----------------------------------------------------------------------*
 |(private)                                                     tk 07/08 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void CONSTRAINTS::MPConstraint2::evaluate_constraint(Teuchos::RCP<Core::FE::Discretization> disc,
    Teuchos::ParameterList& params, Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  if (!(disc->filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!(disc->have_dofs())) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  // see what we have for input
  bool assemblemat1 = systemmatrix1 != Teuchos::null;
  bool assemblemat2 = systemmatrix2 != Teuchos::null;
  bool assemblevec1 = systemvector1 != Teuchos::null;
  bool assemblevec2 = systemvector2 != Teuchos::null;
  bool assemblevec3 = systemvector3 != Teuchos::null;

  // define element matrices and vectors
  Core::LinAlg::SerialDenseMatrix elematrix1;
  Core::LinAlg::SerialDenseMatrix elematrix2;
  Core::LinAlg::SerialDenseVector elevector1;
  Core::LinAlg::SerialDenseVector elevector2;
  Core::LinAlg::SerialDenseVector elevector3;


  const double time = params.get("total time", -1.0);
  const int numcolele = disc->num_my_col_elements();

  // get values from time integrator to scale matrices with
  double scStiff = params.get("scaleStiffEntries", 1.0);
  double scConMat = params.get("scaleConstrMat", 1.0);

  // loop over column elements
  for (int i = 0; i < numcolele; ++i)
  {
    Core::Elements::Element* actele = disc->l_col_element(i);
    Core::Conditions::Condition& cond = *(constrcond_[actele->id()]);
    int condID = cond.parameters().get<int>("ConditionID");

    // computation only if time is larger or equal than initialization time for constraint
    if (inittimes_.find(condID)->second <= time)
    {
      // initialize if it is the first time condition is evaluated
      if (activecons_.find(condID)->second == false)
      {
        const std::string action = params.get<std::string>("action");
        initialize(params, systemvector2);
        params.set("action", action);
      }

      // define global and local index of this bc in redundant vectors
      const int offsetID = params.get<int>("OffsetID");
      int gindex = condID - offsetID;
      const int lindex = (systemvector3->Map()).LID(gindex);

      // Get the current lagrange multiplier value for this condition
      const Teuchos::RCP<Epetra_Vector> lagramul =
          params.get<Teuchos::RCP<Epetra_Vector>>("LagrMultVector");
      const double lagraval = (*lagramul)[lindex];

      // get element location vector, dirichlet flags and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      actele->location_vector(*disc, lm, lmowner, lmstride);
      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();
      if (assemblemat1) elematrix1.shape(eledim, eledim);
      if (assemblemat2) elematrix2.shape(eledim, eledim);
      if (assemblevec1) elevector1.size(eledim);
      if (assemblevec2) elevector2.size(eledim);
      if (assemblevec3) elevector3.size(1);  // elevector3 always contains a scalar

      params.set("ConditionID", condID);
      params.set<Teuchos::RCP<Core::Conditions::Condition>>(
          "condition", Teuchos::rcp(&cond, false));
      // call the element evaluate method
      int err = actele->evaluate(
          params, *disc, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err)
        FOUR_C_THROW(
            "Proc %d: Element %d returned err=%d", disc->get_comm().MyPID(), actele->id(), err);

      int eid = actele->id();

      // Assembly
      if (assemblemat1)
      {
        // scale with time integrator dependent value
        elematrix1.scale(scStiff * lagraval);
        systemmatrix1->assemble(eid, lmstride, elematrix1, lm, lmowner);
      }
      if (assemblemat2)
      {
        std::vector<int> colvec(1);
        colvec[0] = gindex;
        elevector2.scale(scConMat);
        systemmatrix2->assemble(eid, lmstride, elevector2, lm, lmowner, colvec);
      }
      if (assemblevec1)
      {
        elevector1.scale(lagraval);
        Core::LinAlg::Assemble(*systemvector1, elevector1, lm, lmowner);
      }
      if (assemblevec3)
      {
        std::vector<int> constrlm;
        std::vector<int> constrowner;
        constrlm.push_back(gindex);
        constrowner.push_back(actele->owner());
        Core::LinAlg::Assemble(*systemvector3, elevector3, constrlm, constrowner);
      }

      // Load curve business
      const auto* curve = cond.parameters().get_if<int>("curve");
      int curvenum = -1;
      if (curve) curvenum = *curve;
      double curvefac = 1.0;
      bool usetime = true;
      if (time < 0.0) usetime = false;
      if (curvenum >= 0 && usetime)
        curvefac = Global::Problem::instance()
                       ->function_by_id<Core::UTILS::FunctionOfTime>(curvenum)
                       .evaluate(time);
      Teuchos::RCP<Epetra_Vector> timefact =
          params.get<Teuchos::RCP<Epetra_Vector>>("vector curve factors");
      timefact->ReplaceGlobalValues(1, &curvefac, &gindex);
    }
  }
}  // end of evaluate_condition

FOUR_C_NAMESPACE_CLOSE
