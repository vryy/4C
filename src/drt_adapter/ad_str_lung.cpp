/*!----------------------------------------------------------------------

\maintainer Carolin Geitner

\level 3

\brief Structure field adapter for fsi airway simulations with
attached parenchyma balloon

*----------------------------------------------------------------------*/

#include "ad_str_lung.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"
#include "../linalg/linalg_utils.H"

/*======================================================================*/
/* constructor */
ADAPTER::StructureLung::StructureLung(Teuchos::RCP<Structure> stru) : FSIStructureWrapper(stru)
{
  //----------------------------------------------------------------------
  // make sure
  if (structure_ == Teuchos::null) dserror("Failed to create the underlying structural adapter");

  //----------------------------------------------------------------------
  // get lung fluid-structure volume and asi constraints
  std::vector<DRT::Condition*> temp;
  Discretization()->GetCondition("StructFluidSurfCoupling", temp);
  for (unsigned i = 0; i < temp.size(); ++i)
  {
    DRT::Condition& cond = *(temp[i]);
    if (*(cond.Get<std::string>("field")) == "structure")
    {
      constrcond_.push_back(temp[i]);
    }
  }
  if (constrcond_.size() == 0) dserror("No structure-fluid volume constraints found for lung fsi");

  Discretization()->GetCondition("StructAleCoupling", temp);
  for (unsigned i = 0; i < temp.size(); ++i)
  {
    DRT::Condition& cond = *(temp[i]);
    if (*(cond.Get<std::string>("field")) == "structure") asicond_.push_back(temp[i]);
  }
  if (asicond_.size() == 0) dserror("No structure-ale coupling constraints found for lung fsi");

  //----------------------------------------------------------------------
  // consistency check: for each condition, all dofs in the ASI coupling need
  // to be part of the volume coupling, too

  std::map<int, std::set<int>> AllConstrDofMap;
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& constrcond = *(constrcond_[i]);
    int condID = constrcond.GetInt("coupling id");
    const std::vector<int>* constrnodeIDs = constrcond.Nodes();
    std::set<int>& constrdofs = AllConstrDofMap[condID];
    for (unsigned int j = 0; j < constrnodeIDs->size(); ++j)
    {
      int gid = (*constrnodeIDs)[j];
      if (Discretization()->HaveGlobalNode(gid))
      {
        DRT::Node* node = Discretization()->gNode(gid);
        std::vector<int> dofs = Discretization()->Dof(node);
        std::copy(dofs.begin(), dofs.end(), std::inserter(constrdofs, constrdofs.begin()));
      }
    }
  }
  std::map<int, std::set<int>> AllAsiDofMap;
  for (unsigned int i = 0; i < asicond_.size(); ++i)
  {
    DRT::Condition& asicond = *(asicond_[i]);
    int asicondID = asicond.GetInt("coupling id");
    std::set<int>& asidofs = AllAsiDofMap[asicondID];
    const std::vector<int>* asinodeIDs = asicond.Nodes();
    for (unsigned int j = 0; j < asinodeIDs->size(); ++j)
    {
      int gid = (*asinodeIDs)[j];
      if (Discretization()->HaveGlobalNode(gid))
      {
        DRT::Node* node = Discretization()->gNode(gid);
        std::vector<int> dofs = Discretization()->Dof(node);
        std::copy(dofs.begin(), dofs.end(), std::inserter(asidofs, asidofs.begin()));
      }
    }
  }
  for (std::map<int, std::set<int>>::const_iterator it = AllAsiDofMap.begin();
       it != AllAsiDofMap.end(); ++it)
  {
    int condID = it->first;
    const std::set<int>& asidofs = it->second;
    const std::set<int>& constrdofs = AllConstrDofMap[condID];
    std::set<int> intersection;
    std::set_intersection(asidofs.begin(), asidofs.end(), constrdofs.begin(), constrdofs.end(),
        std::inserter(intersection, intersection.begin()));
    if (intersection.size() != asidofs.size())
      dserror(
          "Condition-ID %d: missing ASI dofs or buggy assignment of ASI and volume coupling "
          "condition IDs.",
          condID + 1);
  }

  //----------------------------------------------------------------------
  // build mapextractor for fsi <-> full map
  fsiinterface_ =
      Teuchos::rcp(new LINALG::MapExtractor(*Interface()->FullMap(), Interface()->FSICondMap()));

  //----------------------------------------------------------------------
  // find all dofs belonging to enclosing boundary -> volume coupling dofs
  std::vector<int> dofmapvec;
  std::vector<int> nodes;
  DRT::UTILS::FindConditionedNodes(*Discretization(), "StructFluidSurfCoupling", nodes);
  const int numnode = nodes.size();

  const int ndim = DRT::Problem::Instance()->NDim();
  for (int i = 0; i < numnode; ++i)
  {
    const DRT::Node* actnode = Discretization()->gNode(nodes[i]);
    const std::vector<int> dof = Discretization()->Dof(actnode);
    if (ndim > static_cast<int>(dof.size()))
      dserror(
          "got just %d dofs at node %d (lid=%d) but expected %d", dof.size(), nodes[i], i, ndim);
    std::copy(&dof[0], &dof[0] + ndim, back_inserter(dofmapvec));
  }

  lungconstraintmap_ = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0], 0, Discretization()->Comm()));
}


/*======================================================================*/
/* list of coupled fluid-structure volumes */
void ADAPTER::StructureLung::ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID)
{
  MinLungVolConID = 1;

  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);
    int condID = cond.GetInt("coupling id");
    if (LungVolConIDs.find(condID) == LungVolConIDs.end())
    {
      if (condID < MinLungVolConID) MinLungVolConID = condID;
      LungVolConIDs.insert(condID);
    }
  }
}


/*======================================================================*/
/* determine initial volumes */
void ADAPTER::StructureLung::InitializeVolCon(
    Teuchos::RCP<Epetra_Vector> initvol, Teuchos::RCP<Epetra_Vector> signvol, const int offsetID)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  if (initvol == Teuchos::null or signvol == Teuchos::null)
    dserror("Structure lung volume constraint cannot be initialized");

  if (!initvol->Map().SameAs(signvol->Map()))
    dserror("Structure lung volume constraint cannot be initialized");

  Teuchos::ParameterList params;
  params.set("action", "calc_struct_constrvol");

  // set displacements
  Discretization()->ClearState();
  Discretization()->SetState("displacement", Dispnp());

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    {
      // Get ConditionID of current condition if defined and write value in parameterlist
      int condID = cond.GetInt("coupling id");
      params.set("ConditionID", condID);
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
        curr->second->LocationVector(*Discretization(), lm, lmowner, lmstride);

        // reshape element matrices and vectors and init to zero
        elevector3.Size(1);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params, *Discretization(), lm, elematrix1, elematrix2,
            elevector1, elevector2, elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly

        std::vector<int> constrlm;
        std::vector<int> constrowner;
        constrlm.push_back(condID - offsetID);
        constrowner.push_back(curr->second->Owner());
        LINALG::Assemble(*initvol, elevector3, constrlm, constrowner);
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
void ADAPTER::StructureLung::EvaluateVolCon(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> StructMatrix, Teuchos::RCP<Epetra_Vector> StructRHS,
    Teuchos::RCP<Epetra_Vector> CurrVols, Teuchos::RCP<Epetra_Vector> SignVols,
    Teuchos::RCP<Epetra_Vector> lagrMultVecRed, const int offsetID)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // parameter list
  Teuchos::ParameterList params;
  params.set("action", "calc_struct_volconstrstiff");

  // set displacements
  Discretization()->ClearState();
  Discretization()->SetState("displacement", Dispnp());

  //---------------------------------------------------------------------
  // loop through conditions and evaluate them
  //---------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond.GetInt("coupling id");
    params.set("ConditionID", condID);

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition>>("condition", Teuchos::rcp(&cond, false));

    // global and local ID of this bc in the redundant vectors
    const int gindex = condID - offsetID;
    const int lindex = (CurrVols->Map()).LID(gindex);
    if (lindex == -1) dserror("Corrupt vector of current volumes");

    // Get the current lagrange multiplier value for this condition
    const double lagraval = (*lagrMultVecRed)[lindex];

    // Get the sign (of volume) for this condition
    const double sign = (*SignVols)[lindex];

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // there might be processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;

    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*Discretization(), lm, lmowner, lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();
      elematrix1.Shape(eledim, eledim);  // stiffness part
      elematrix2.Shape(
          eledim, eledim);      // this element matrix is only needed for the function call only
      elevector1.Size(eledim);  // rhs part
      elevector2.Size(eledim);  // constraint matrix
      elevector3.Size(1);       // current volume

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params, *Discretization(), lm, elematrix1, elematrix2,
          elevector1, elevector2, elevector3);
      if (err) dserror("error while evaluating elements");

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

      elematrix1.Scale(-lagraval * sign);
      StructMatrix->Assemble(eid, lmstride, elematrix1, lm, lmowner);

      // assemble to rectangular matrix. The column corresponds to the constraint ID.
      std::vector<int> colvec(1);
      colvec[0] = gindex;
      elevector2.Scale(-sign);
      StructMatrix->Assemble(eid, lmstride, elevector2, lm, lmowner, colvec);

      // "Newton-ready" residual -> already scaled with -1.0
      elevector1.Scale(lagraval * sign);
      LINALG::Assemble(*StructRHS, elevector1, lm, lmowner);

      // No scaling with -1.0 necessary here, since the constraint rhs is determined consistently,
      // i.e.  -(Vcurr - Vold) in the fsi algorithm, thus -1.0 is included there.
      std::vector<int> constrlm;
      std::vector<int> constrowner;
      constrlm.push_back(gindex);
      constrowner.push_back(curr->second->Owner());
      elevector3.Scale(sign);
      LINALG::Assemble(*CurrVols, elevector3, constrlm, constrowner);
    }
  }

  StructMatrix->Complete();
  StructMatrix->Matrix(1, 0) = *StructMatrix->Matrix(0, 1).Transpose();

#if 0   // Debug
  //############################################################################################
  // test
  LINALG::SparseMatrix& StructConstrMatrix = StructMatrix->Matrix(0,1);
  Teuchos::RCP<Epetra_CrsMatrix> Aprime = StructConstrMatrix.EpetraMatrix();
  int MaxNumEntries = Aprime->MaxNumEntries();
  int NumEntries;

  std::set<int> ConIDs;
  for (unsigned int i=0; i<constrcond_.size();++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);
    int condID = cond.GetInt("coupling id");
    if (ConIDs.find(condID) == ConIDs.end())
    {
      ConIDs.insert(condID);
    }
  }
  std::map<int,std::vector<int> > AllColIndices;
  for (std::set<int>::const_iterator it = ConIDs.begin(); it != ConIDs.end(); ++it)
  {
    // check StructConstrMat
    int condID = *it;
    const int gindex = condID-offsetID;
    std::vector<int> & ColIndices = AllColIndices[condID];

    for (int j=0; j<Aprime->NumMyRows(); ++j)
    {
      std::vector<int>    Indices(MaxNumEntries);
      std::vector<double> Values(MaxNumEntries);
      int Row = Aprime->GRID(j);
      int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
      if (ierr) dserror("ExtractGlobalRowCopy failed: err=%d",ierr);
      for (unsigned int k=0; k<NumEntries; ++k)
      {
        if (Indices[k] == gindex)
          ColIndices.push_back(Row);
      }
    }
  }

  std::map<int,std::vector<int> > AllAsiDofMap;
  for (unsigned int i = 0; i < asicond_.size(); ++i)
  {
    DRT::Condition& asicond = *(asicond_[i]);
    int asicondID = asicond.GetInt("coupling id");
    std::vector<int> &  asidofs  = AllAsiDofMap[asicondID];
    const std::vector<int>* asinodeIDs = asicond.Nodes();
    for (unsigned int j = 0; j < asinodeIDs->size(); ++j)
    {
      int gid = (*asinodeIDs)[j];
      if (Discretization()->HaveGlobalNode(gid))
      {
        DRT::Node* node = Discretization()->gNode(gid);
        if (node->Owner() == Discretization()->Comm().MyPID())
        {
          std::vector<int> dofs = Discretization()->Dof(node);
          for (unsigned int k = 0; k < dofs.size(); ++k)
            asidofs.push_back(dofs[k]);
        }
      }
    }
  }

  for (std::map<int,std::vector<int> >::const_iterator it=AllAsiDofMap.begin(); it!=AllAsiDofMap.end(); ++it)
  {
    int condID = it->first;
    const std::vector<int> & asidofs = it->second;
    const std::vector<int> & colindices = AllColIndices[condID];

    for (unsigned int i=0; i<asidofs.size(); ++i)
    {
      bool found = false;
      for (unsigned int j=0; j<colindices.size(); ++j)
      {
        if (asidofs[i] == colindices[j])
        {
          found = true;
          break;
        }
      }
      if (found == false)
      {
        dserror("something strange in StructConstrMatrix: missing asi dof %d",asidofs[i]);
      }
    }
  }

  // check ConstrStructMat

  for (std::map<int,std::vector<int> >::const_iterator it=AllAsiDofMap.begin(); it!=AllAsiDofMap.end(); ++it)
  {
    int condID = it->first;
    std::vector<int> asidofs = it->second;
    std::vector<int> colindices = AllColIndices[condID];
    int numcolind = colindices.size();
    Epetra_Map ColIndices(-1, numcolind, &colindices[0], 0, Discretization()->Comm());
    Epetra_Map RedColIndices = *(LINALG::AllreduceEMap(ColIndices));

    for (unsigned int k=0; k<asidofs.size(); ++k)
    {
      int myind = asidofs[k];
      if (RedColIndices.LID(myind) == -1)
        dserror("ColIndices and Indices have different entries!");
    }
  }
  // end test
  //############################################################################################
#endif  // End debug


  // Apply Dirichlet BC to stiffness matrix and rhs vector and
  // exclude all forces on the outflow boundary (except those at the fsi
  // partition) and corresponding stiffness matrix contributions
  const Teuchos::RCP<const Epetra_Map>& condmap = GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map>& outflowmap = Interface()->Map(2);
  Teuchos::RCP<Epetra_Map> finmap = LINALG::MergeMap(*condmap, *outflowmap, false);
  StructMatrix->ApplyDirichlet(*finmap, false);

  const Teuchos::RCP<const Epetra_Map>& dispmap = StructMatrix->RangeExtractor().Map(0);
  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dispmap, true);
  GetDBCMapExtractor()->InsertCondVector(GetDBCMapExtractor()->ExtractCondVector(zeros), StructRHS);
  Interface()->InsertVector(Interface()->ExtractVector(zeros, 2), 2, StructRHS);
}


/*======================================================================*/
/* write additional stuff in case of restart */
void ADAPTER::StructureLung::WriteVolConRestart(Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
    Teuchos::RCP<Epetra_Vector> OldVolsRed, Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
{
  Teuchos::RCP<IO::DiscretizationWriter> output = DiscWriter();
  std::stringstream stream1;
  stream1 << "OldVols";
  std::stringstream stream2;
  stream2 << "OldFlowRates";
  std::stringstream stream3;
  stream3 << "OldLagrMult";

  // NOTE: OldFlowRatesRed, OldVolsRed and OldLagrMultRed are redundant vectors,
  // hence "conversion" into a std::vector<double> provides identical
  // results on all processors. However, only processor 0 writes
  // output in WriteRedundantDoubleVector.

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

  output->WriteRedundantDoubleVector(stream1.str(), volumes);
  output->WriteRedundantDoubleVector(stream2.str(), flowrates);
  output->WriteRedundantDoubleVector(stream3.str(), lmult);
}


/*======================================================================*/
/* output of volume constraint related forces*/
void ADAPTER::StructureLung::OutputForces(Teuchos::RCP<Epetra_Vector> Forces)
{
  Teuchos::RCP<IO::DiscretizationWriter> output = DiscWriter();
  output->WriteVector("Add_Forces", Forces);
}


/*======================================================================*/
/* read additional stuff in case of restart */
void ADAPTER::StructureLung::ReadVolConRestart(const int step,
    Teuchos::RCP<Epetra_Vector> OldFlowRatesRed, Teuchos::RCP<Epetra_Vector> OldVolsRed,
    Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
{
  IO::DiscretizationReader reader(Discretization(), step);
  if (step != reader.ReadInt("step")) dserror("Time step on file not equal to given step");
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

  reader.ReadRedundantDoubleVector(volumes, stream1.str());
  reader.ReadRedundantDoubleVector(flowrates, stream2.str());
  reader.ReadRedundantDoubleVector(lmult, stream3.str());

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
