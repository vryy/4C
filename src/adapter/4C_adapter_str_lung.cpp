/*----------------------------------------------------------------------*/
/*! \file


\level 3

\brief Structure field adapter for fsi airway simulations with
attached parenchyma balloon

*----------------------------------------------------------------------*/

#include "4C_adapter_str_lung.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::StructureLung::StructureLung(Teuchos::RCP<Structure> stru) : FSIStructureWrapper(stru)
{
  //----------------------------------------------------------------------
  // make sure
  if (structure_ == Teuchos::null)
    FOUR_C_THROW("Failed to create the underlying structural adapter");

  //----------------------------------------------------------------------
  // get lung fluid-structure volume and asi constraints
  std::vector<Core::Conditions::Condition*> temp;
  discretization()->GetCondition("StructFluidSurfCoupling", temp);
  for (unsigned i = 0; i < temp.size(); ++i)
  {
    Core::Conditions::Condition& cond = *(temp[i]);
    if ((cond.parameters().get<std::string>("field")) == "structure")
    {
      constrcond_.push_back(temp[i]);
    }
  }
  if (constrcond_.size() == 0)
    FOUR_C_THROW("No structure-fluid volume constraints found for lung fsi");

  discretization()->GetCondition("StructAleCoupling", temp);
  for (unsigned i = 0; i < temp.size(); ++i)
  {
    Core::Conditions::Condition& cond = *(temp[i]);
    if ((cond.parameters().get<std::string>("field")) == "structure") asicond_.push_back(temp[i]);
  }
  if (asicond_.size() == 0)
    FOUR_C_THROW("No structure-ale coupling constraints found for lung fsi");

  //----------------------------------------------------------------------
  // consistency check: for each condition, all dofs in the ASI coupling need
  // to be part of the volume coupling, too

  std::map<int, std::set<int>> AllConstrDofMap;
  for (auto& constrcond : constrcond_)
  {
    int condID = constrcond->parameters().get<int>("coupling id");
    const std::vector<int>* constrnodeIDs = constrcond->GetNodes();
    std::set<int>& constrdofs = AllConstrDofMap[condID];
    for (int gid : *constrnodeIDs)
    {
      if (discretization()->HaveGlobalNode(gid))
      {
        Core::Nodes::Node* node = discretization()->gNode(gid);
        std::vector<int> dofs = discretization()->Dof(node);
        std::copy(dofs.begin(), dofs.end(), std::inserter(constrdofs, constrdofs.begin()));
      }
    }
  }
  std::map<int, std::set<int>> AllAsiDofMap;
  for (auto& asicond : asicond_)
  {
    int asicondID = asicond->parameters().get<int>("coupling id");
    std::set<int>& asidofs = AllAsiDofMap[asicondID];
    const std::vector<int>* asinodeIDs = asicond->GetNodes();
    for (int gid : *asinodeIDs)
    {
      if (discretization()->HaveGlobalNode(gid))
      {
        Core::Nodes::Node* node = discretization()->gNode(gid);
        std::vector<int> dofs = discretization()->Dof(node);
        std::copy(dofs.begin(), dofs.end(), std::inserter(asidofs, asidofs.begin()));
      }
    }
  }
  for (const auto& [condID, asidofs] : AllAsiDofMap)
  {
    const std::set<int>& constrdofs = AllConstrDofMap[condID];
    std::set<int> intersection;
    std::set_intersection(asidofs.begin(), asidofs.end(), constrdofs.begin(), constrdofs.end(),
        std::inserter(intersection, intersection.begin()));
    if (intersection.size() != asidofs.size())
      FOUR_C_THROW(
          "Condition-ID %d: missing ASI dofs or buggy assignment of ASI and volume coupling "
          "condition IDs.",
          condID + 1);
  }

  //----------------------------------------------------------------------
  // build mapextractor for fsi <-> full map
  fsiinterface_ = Teuchos::rcp(
      new Core::LinAlg::MapExtractor(*Interface()->FullMap(), Interface()->FSICondMap()));

  //----------------------------------------------------------------------
  // find all dofs belonging to enclosing boundary -> volume coupling dofs
  std::vector<int> dofmapvec;
  std::vector<int> nodes;
  Core::Conditions::FindConditionedNodes(*discretization(), "StructFluidSurfCoupling", nodes);
  const int numnode = nodes.size();

  const int ndim = Global::Problem::Instance()->NDim();
  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* actnode = discretization()->gNode(nodes[i]);
    const std::vector<int> dof = discretization()->Dof(actnode);
    if (ndim > static_cast<int>(dof.size()))
      FOUR_C_THROW(
          "got just %d dofs at node %d (lid=%d) but expected %d", dof.size(), nodes[i], i, ndim);
    std::copy(dof.data(), dof.data() + ndim, back_inserter(dofmapvec));
  }

  lungconstraintmap_ = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec.size(), dofmapvec.data(), 0, discretization()->Comm()));
}


/*======================================================================*/
/* list of coupled fluid-structure volumes */
void Adapter::StructureLung::ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID)
{
  MinLungVolConID = 1;

  for (auto& cond : constrcond_)
  {
    int condID = cond->parameters().get<int>("coupling id");
    if (LungVolConIDs.find(condID) == LungVolConIDs.end())
    {
      if (condID < MinLungVolConID) MinLungVolConID = condID;
      LungVolConIDs.insert(condID);
    }
  }
}


/*======================================================================*/
/* determine initial volumes */
void Adapter::StructureLung::InitializeVolCon(
    Teuchos::RCP<Epetra_Vector> initvol, Teuchos::RCP<Epetra_Vector> signvol, const int offsetID)
{
  if (!(discretization()->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!discretization()->HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  if (initvol == Teuchos::null or signvol == Teuchos::null)
    FOUR_C_THROW("Structure lung volume constraint cannot be initialized");

  if (!initvol->Map().SameAs(signvol->Map()))
    FOUR_C_THROW("Structure lung volume constraint cannot be initialized");

  Teuchos::ParameterList params;
  params.set("action", "calc_struct_constrvol");

  // set displacements
  discretization()->ClearState();
  discretization()->set_state("displacement", Dispnp());

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    Core::Conditions::Condition& cond = *(constrcond_[i]);

    {
      // Get ConditionID of current condition if defined and write value in parameterlist
      int condID = cond.parameters().get<int>("coupling id");
      params.set("ConditionID", condID);
      params.set<Teuchos::RCP<Core::Conditions::Condition>>(
          "condition", Teuchos::rcp(&cond, false));

      // define element matrices and vectors
      Core::LinAlg::SerialDenseMatrix elematrix1;
      Core::LinAlg::SerialDenseMatrix elematrix2;
      Core::LinAlg::SerialDenseVector elevector1;
      Core::LinAlg::SerialDenseVector elevector2;
      Core::LinAlg::SerialDenseVector elevector3;

      std::map<int, Teuchos::RCP<Core::Elements::Element>>& geom = cond.Geometry();
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
        curr->second->LocationVector(*discretization(), lm, lmowner, lmstride);

        // reshape element matrices and vectors and init to zero
        elevector3.size(1);

        // call the element specific evaluate method
        int err = curr->second->evaluate(params, *discretization(), lm, elematrix1, elematrix2,
            elevector1, elevector2, elevector3);
        if (err) FOUR_C_THROW("error while evaluating elements");

        // assembly

        std::vector<int> constrlm;
        std::vector<int> constrowner;
        constrlm.push_back(condID - offsetID);
        constrowner.push_back(curr->second->Owner());
        Core::LinAlg::Assemble(*initvol, elevector3, constrlm, constrowner);
      }
    }
  }

  // determine sign of all calculated volumes in _redundant_ vector
  //
  // the functions determining current volumes may provide also negative
  // volumes. this is ok for a "normal" volume constraint where the
  // difference between a reference volume and the current one (both
  // determined the same way, thus with the same sign) is
  // needed. since this difference is required to be zero anyway, the
  // sign does not play a role. however, in the lung volume
  // constraint, the change in structural volume needs to be
  // correlated with the outflowing fluid volume. hence, the sign of
  // this change (volume decrease or increase!) _does_
  // matter. therefore, we check whether the initial structural volume
  // calculated for each constraint is positive or negative, save this
  // sign and multiply all structural volume constraint related stuff
  // with this sign later on.

  for (int i = 0; i < initvol->MyLength(); ++i)
  {
    double sumvol = 0.;
    signvol->Map().Comm().SumAll(&((*initvol)[i]), &sumvol, 1);
    if (sumvol > 0.0)
    {
      (*signvol)[i] = 1.0;
    }
    else
    {
      (*signvol)[i] = -1.0;
      (*initvol)[i] *= -1.0;
    }
  }
}


/*======================================================================*/
/* evaluate structural part of fluid-structure volume constraint */
void Adapter::StructureLung::EvaluateVolCon(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> StructMatrix,
    Teuchos::RCP<Epetra_Vector> StructRHS, Teuchos::RCP<Epetra_Vector> CurrVols,
    Teuchos::RCP<Epetra_Vector> SignVols, Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
    const int offsetID)
{
  if (!(discretization()->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!discretization()->HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  // parameter list
  Teuchos::ParameterList params;
  params.set("action", "calc_struct_volconstrstiff");

  // set displacements
  discretization()->ClearState();
  discretization()->set_state("displacement", Dispnp());

  //---------------------------------------------------------------------
  // loop through conditions and evaluate them
  //---------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    Core::Conditions::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond.parameters().get<int>("coupling id");
    params.set("ConditionID", condID);

    // elements might need condition
    params.set<Teuchos::RCP<Core::Conditions::Condition>>("condition", Teuchos::rcp(&cond, false));

    // global and local ID of this bc in the redundant vectors
    const int gindex = condID - offsetID;
    const int lindex = (CurrVols->Map()).LID(gindex);
    if (lindex == -1) FOUR_C_THROW("Corrupt vector of current volumes");

    // Get the current lagrange multiplier value for this condition
    const double lagraval = (*lagrMultVecRed)[lindex];

    // Get the sign (of volume) for this condition
    const double sign = (*SignVols)[lindex];

    // define element matrices and vectors
    Core::LinAlg::SerialDenseMatrix elematrix1;
    Core::LinAlg::SerialDenseMatrix elematrix2;
    Core::LinAlg::SerialDenseVector elevector1;
    Core::LinAlg::SerialDenseVector elevector2;
    Core::LinAlg::SerialDenseVector elevector3;

    std::map<int, Teuchos::RCP<Core::Elements::Element>>& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // there might be processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator curr;

    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*discretization(), lm, lmowner, lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();
      elematrix1.shape(eledim, eledim);  // stiffness part
      elematrix2.shape(
          eledim, eledim);      // this element matrix is only needed for the function call only
      elevector1.size(eledim);  // rhs part
      elevector2.size(eledim);  // constraint matrix
      elevector3.size(1);       // current volume

      // call the element specific evaluate method
      int err = curr->second->evaluate(params, *discretization(), lm, elematrix1, elematrix2,
          elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("error while evaluating elements");

      // assembly
      int eid = curr->second->Id();

      // distinction whether this part of the condition belongs to the inflow boundary (in this case
      // only contribution to overall volume and corresponding constraint matrix) or to the rest (in
      // this case additional rhs and matrix contributions in rows corresponding to structural dofs)

      // NOTE:
      // - no time integration related scaling needed here, since everything is evaluated at the
      //   end of the time step and assembled to the already completed system.
      // - all element quantities are multiplied by the previously derived sign in order to assure
      //   that computed volumes are positive and corresponding derivatives are determined with a
      //   consistent sign.
      // - additional multiplication with -1.0 accounts for different definition of volume
      //   difference ("normal" volume constraint: Vref-Vcurr, lung volume constraint: Vcurr-Vref,
      //   which seems more natural in this case)

      elematrix1.scale(-lagraval * sign);
      StructMatrix->Assemble(eid, lmstride, elematrix1, lm, lmowner);

      // assemble to rectangular matrix. The column corresponds to the constraint ID.
      std::vector<int> colvec(1);
      colvec[0] = gindex;
      elevector2.scale(-sign);
      StructMatrix->Assemble(eid, lmstride, elevector2, lm, lmowner, colvec);

      // "Newton-ready" residual -> already scaled with -1.0
      elevector1.scale(lagraval * sign);
      Core::LinAlg::Assemble(*StructRHS, elevector1, lm, lmowner);

      // No scaling with -1.0 necessary here, since the constraint rhs is determined consistently,
      // i.e.  -(Vcurr - Vold) in the fsi algorithm, thus -1.0 is included there.
      std::vector<int> constrlm;
      std::vector<int> constrowner;
      constrlm.push_back(gindex);
      constrowner.push_back(curr->second->Owner());
      elevector3.scale(sign);
      Core::LinAlg::Assemble(*CurrVols, elevector3, constrlm, constrowner);
    }
  }

  StructMatrix->Complete();
  StructMatrix->Matrix(1, 0) = *StructMatrix->Matrix(0, 1).Transpose();

  // Apply Dirichlet BC to stiffness matrix and rhs vector and
  // exclude all forces on the outflow boundary (except those at the fsi
  // partition) and corresponding stiffness matrix contributions
  const Teuchos::RCP<const Epetra_Map>& condmap = GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map>& outflowmap = Interface()->Map(2);
  Teuchos::RCP<Epetra_Map> finmap = Core::LinAlg::MergeMap(*condmap, *outflowmap, false);
  StructMatrix->ApplyDirichlet(*finmap, false);

  const Teuchos::RCP<const Epetra_Map>& dispmap = StructMatrix->RangeExtractor().Map(0);
  Teuchos::RCP<Epetra_Vector> zeros = Core::LinAlg::CreateVector(*dispmap, true);
  GetDBCMapExtractor()->InsertCondVector(GetDBCMapExtractor()->ExtractCondVector(zeros), StructRHS);
  Interface()->InsertVector(Interface()->ExtractVector(zeros, 2), 2, StructRHS);
}


/*======================================================================*/
/* write additional stuff in case of restart */
void Adapter::StructureLung::WriteVolConRestart(Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
    Teuchos::RCP<Epetra_Vector> OldVolsRed, Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
{
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = disc_writer();
  std::stringstream stream1;
  stream1 << "OldVols";
  std::stringstream stream2;
  stream2 << "OldFlowRates";
  std::stringstream stream3;
  stream3 << "OldLagrMult";

  // NOTE: OldFlowRatesRed, OldVolsRed and OldLagrMultRed are redundant vectors,
  // hence "conversion" into a std::vector<double> provides identical
  // results on all processors. However, only processor 0 writes
  // output in write_redundant_double_vector.

  Teuchos::RCP<std::vector<double>> flowrates =
      Teuchos::rcp(new std::vector<double>(OldFlowRatesRed->MyLength()));
  for (int i = 0; i < OldFlowRatesRed->MyLength(); ++i)
  {
    (*flowrates)[i] = (*OldFlowRatesRed)[i];
  }
  Teuchos::RCP<std::vector<double>> volumes =
      Teuchos::rcp(new std::vector<double>(OldVolsRed->MyLength()));
  for (int i = 0; i < OldVolsRed->MyLength(); ++i)
  {
    (*volumes)[i] = (*OldVolsRed)[i];
  }
  Teuchos::RCP<std::vector<double>> lmult =
      Teuchos::rcp(new std::vector<double>(OldLagrMultRed->MyLength()));
  for (int i = 0; i < OldLagrMultRed->MyLength(); ++i)
  {
    (*lmult)[i] = (*OldLagrMultRed)[i];
  }

  output->write_redundant_double_vector(stream1.str(), volumes);
  output->write_redundant_double_vector(stream2.str(), flowrates);
  output->write_redundant_double_vector(stream3.str(), lmult);
}


/*======================================================================*/
/* output of volume constraint related forces*/
void Adapter::StructureLung::OutputForces(Teuchos::RCP<Epetra_Vector> Forces)
{
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = disc_writer();
  output->write_vector("Add_Forces", Forces);
}


/*======================================================================*/
/* read additional stuff in case of restart */
void Adapter::StructureLung::ReadVolConRestart(const int step,
    Teuchos::RCP<Epetra_Vector> OldFlowRatesRed, Teuchos::RCP<Epetra_Vector> OldVolsRed,
    Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
{
  Core::IO::DiscretizationReader reader(
      discretization(), Global::Problem::Instance()->InputControlFile(), step);
  if (step != reader.read_int("step")) FOUR_C_THROW("Time step on file not equal to given step");
  std::stringstream stream1;
  stream1 << "OldVols";
  std::stringstream stream2;
  stream2 << "OldFlowRates";
  std::stringstream stream3;
  stream3 << "OldLagrMult";

  Teuchos::RCP<std::vector<double>> flowrates =
      Teuchos::rcp(new std::vector<double>(OldFlowRatesRed->MyLength()));
  Teuchos::RCP<std::vector<double>> volumes =
      Teuchos::rcp(new std::vector<double>(OldVolsRed->MyLength()));
  Teuchos::RCP<std::vector<double>> lmult =
      Teuchos::rcp(new std::vector<double>(OldLagrMultRed->MyLength()));

  reader.read_redundant_double_vector(volumes, stream1.str());
  reader.read_redundant_double_vector(flowrates, stream2.str());
  reader.read_redundant_double_vector(lmult, stream3.str());

  for (int i = 0; i < OldFlowRatesRed->MyLength(); ++i)
  {
    (*OldFlowRatesRed)[i] = (*flowrates)[i];
  }
  for (int i = 0; i < OldVolsRed->MyLength(); ++i)
  {
    (*OldVolsRed)[i] = (*volumes)[i];
  }
  for (int i = 0; i < OldLagrMultRed->MyLength(); ++i)
  {
    (*OldLagrMultRed)[i] = (*lmult)[i];
  }
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
