/*!----------------------------------------------------------------------

\brief Fluid field adapter for fsi airway simulations with attached
parenchyma balloon

\maintainer Carolin Geitner



*----------------------------------------------------------------------*/

#include "ad_fld_lung.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_io/io.H"
#include "../drt_fluid_ele/fluid_ele_action.H"


/*======================================================================*/
/* constructor */
ADAPTER::FluidLung::FluidLung(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure
  if (fluid_ == Teuchos::null) dserror("Failed to create the underlying fluid adapter");
  return;
}

/*======================================================================*/
/* initialization */
void ADAPTER::FluidLung::Init()
{
  // call base class init
  FluidFSI::Init();

  // get lung fluid-structure volume constraints

  std::vector<DRT::Condition*> temp;
  Discretization()->GetCondition("StructFluidSurfCoupling", temp);
  for (unsigned i = 0; i < temp.size(); ++i)
  {
    DRT::Condition& cond = *(temp[i]);
    if (*(cond.Get<std::string>("field")) == "fluid") constrcond_.push_back(temp[i]);
  }
  if (constrcond_.size() == 0) dserror("No structure-fluid volume constraints found for lung fsi");

  // build map extractor for fsi <-> full map

  fsiinterface_ =
      Teuchos::rcp(new LINALG::MapExtractor(*Interface()->FullMap(), Interface()->FSICondMap()));

  // build map extractor for asi, other <-> full inner map

  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(Interface()->OtherMap());
  maps.push_back(Interface()->LungASICondMap());
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  innersplit_ = Teuchos::rcp(new LINALG::MapExtractor(*fullmap, Interface()->LungASICondMap()));

  // build mapextractor for outflow fsi boundary dofs <-> full map

  std::vector<int> fsinodes;
  DRT::UTILS::FindConditionedNodes(*Discretization(), "FSICoupling", fsinodes);

  std::set<int> outflownodes;
  DRT::UTILS::FindConditionedNodes(*Discretization(), "StructAleCoupling", outflownodes);

  std::vector<int> outflowfsinodes;

  for (unsigned int i = 0; i < fsinodes.size(); ++i)
  {
    if (outflownodes.find(fsinodes[i]) != outflownodes.end())
      outflowfsinodes.push_back(fsinodes[i]);
  }

  std::vector<int> dofmapvec;

  for (unsigned int i = 0; i < outflowfsinodes.size(); ++i)
  {
    DRT::Node* actnode = Discretization()->gNode(outflowfsinodes[i]);
    const std::vector<int> dof = Discretization()->Dof(actnode);

    const int ndim = DRT::Problem::Instance()->NDim();
    if (ndim > static_cast<int>(dof.size()))
      dserror("got just %d dofs but expected %d", dof.size(), ndim);
    std::copy(&dof[0], &dof[0] + ndim, back_inserter(dofmapvec));
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(), dofmapvec.end());
  if (pos != dofmapvec.end() and *pos < 0) dserror("illegal dof number %d", *pos);

  Teuchos::RCP<Epetra_Map> outflowfsidofmap = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0], 0, Discretization()->Comm()));

  outflowfsiinterface_ =
      Teuchos::rcp(new LINALG::MapExtractor(*Interface()->FullMap(), outflowfsidofmap));
}


/*======================================================================*/
/* list of coupled fluid-structure volumes */
void ADAPTER::FluidLung::ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID)
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
/* determine initial flow rates */
void ADAPTER::FluidLung::InitializeVolCon(
    Teuchos::RCP<Epetra_Vector> initflowrate, const int offsetID)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // set ale displacements, fluid and grid velocities
  Discretization()->ClearState();
  Discretization()->SetState("convectivevel", ConvectiveVel());
  Discretization()->SetState("dispnp", Dispnp());

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist

    int condID = cond.GetInt("coupling id");

    Teuchos::ParameterList params;
    params.set("ConditionID", condID);
    params.set<Teuchos::RCP<DRT::Condition>>("condition", Teuchos::rcp(&cond, false));
    params.set<int>("action", FLD::flowratederiv);
    params.set("flowrateonly", true);
    const double dt = Dt();
    params.set("dt", dt);

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

      // Reshape element matrices and vectors and init to zero
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
      LINALG::Assemble(*initflowrate, elevector3, constrlm, constrowner);
    }
  }
}


/*======================================================================*/
/* evaluate structural part of fluid-structure volume constraint */
void ADAPTER::FluidLung::EvaluateVolCon(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> FluidShapeDerivMatrix,
    Teuchos::RCP<LINALG::SparseMatrix> FluidConstrMatrix,
    Teuchos::RCP<LINALG::SparseMatrix> ConstrFluidMatrix,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> AleConstrMatrix,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ConstrAleMatrix,
    Teuchos::RCP<Epetra_Vector> FluidRHS, Teuchos::RCP<Epetra_Vector> CurrFlowRates,
    Teuchos::RCP<Epetra_Vector> lagrMultVecRed, const int offsetID, const double dttheta)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // set ale displacements, fluid and grid velocities
  Discretization()->ClearState();
  Discretization()->SetState("convectivevel", ConvectiveVel());
  Discretization()->SetState("dispnp", Dispnp());

  // residual scaling for fluid matrices needed
  const double invresscale = 1.0 / ResidualScaling();

  //---------------------------------------------------------------------
  // loop through conditions and evaluate them
  //---------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond.GetInt("coupling id");
    Teuchos::ParameterList params;
    params.set("ConditionID", condID);
    const double dt = Dt();
    params.set("dt", dt);

    // global and local ID of this bc in the redundant vectors
    int gindex = condID - offsetID;
    const int lindex = (CurrFlowRates->Map()).LID(gindex);

    // Get the current lagrange multiplier value for this condition
    const double lagraval = (*lagrMultVecRed)[lindex];

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition>>("condition", Teuchos::rcp(&cond, false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;  // (d^2 Q)/(du dd)
    Epetra_SerialDenseMatrix elematrix2;  // (d^2 Q)/(dd)^2
    Epetra_SerialDenseVector elevector1;  // dQ/du
    Epetra_SerialDenseVector elevector2;  // dQ/dd
    Epetra_SerialDenseVector elevector3;  // Q

    std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // there might be processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;

    // define element action
    params.set<int>("action", FLD::flowratederiv);

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
      elematrix1.Shape(eledim, eledim);
      elematrix2.Shape(eledim, eledim);
      elevector1.Size(eledim);
      elevector2.Size(eledim);
      elevector3.Size(1);

      //---------------------------------------------------------------------
      // call the element specific evaluate method
      int err = curr->second->Evaluate(params, *Discretization(), lm, elematrix1, elematrix2,
          elevector1, elevector2, elevector3);
      if (err) dserror("error while evaluating elements");

      //---------------------------------------------------------------------
      // assembly
      int eid = curr->second->Id();

      elematrix1.Scale(-lagraval * invresscale);
      FluidShapeDerivMatrix->Assemble(eid, lmstride, elematrix1, lm, lmowner);

      // assemble to rectangular matrix. The column corresponds to the constraint ID.
      std::vector<int> colvec(1);
      colvec[0] = gindex;

      elevector1.Scale(-1.0);
      FluidConstrMatrix->Assemble(eid, elevector1, lm, lmowner, colvec);

      elevector2.Scale(-dttheta);
      AleConstrMatrix->Assemble(eid, lmstride, elevector2, lm, lmowner, colvec);

      // negative sign (for shift to rhs) is already implicitly taken into account!
      elevector1.Scale(-lagraval * invresscale);
      LINALG::Assemble(*FluidRHS, elevector1, lm, lmowner);

      std::vector<int> constrlm;
      std::vector<int> constrowner;
      constrlm.push_back(gindex);
      constrowner.push_back(curr->second->Owner());
      LINALG::Assemble(*CurrFlowRates, elevector3, constrlm, constrowner);
    }
  }

  FluidShapeDerivMatrix->Complete();

  const Epetra_Map& constrmap = ConstrFluidMatrix->RangeMap();
  const Epetra_Map& fluidmap = FluidConstrMatrix->RangeMap();
  FluidConstrMatrix->Complete(constrmap, fluidmap);
  AleConstrMatrix->Complete();

  // transposed "ale" constraint matrix -> linearization of constraint equation
  for (int i = 0; i < Interface()->NumMaps(); ++i)
    ConstrAleMatrix->Matrix(0, i) = *AleConstrMatrix->Matrix(i, 0).Transpose();
  ConstrAleMatrix->Complete();

  // Note: there is no contribution of the FSI boundary to the overall
  // outflow, since u-u_grid is always zero here. However, in the above
  // calculation of corresponding matrices and vectors, the coupling of
  // displacements and velocities at the interface was not considered (i.e. a
  // variation of the flowrate w.r.t. the velocities always involves also a
  // variation w.r.t. the dependent displacements). Therefore, the
  // corresponding contributions are not zero here and need to be set to zero
  // in the following.

  // At the outlet, no fluid Dirichlet conditions are present!

  const Teuchos::RCP<const Epetra_Map>& outflowfsimap = outflowfsiinterface_->Map(1);
  FluidShapeDerivMatrix->ApplyDirichlet(*outflowfsimap, false);
  FluidConstrMatrix->ApplyDirichlet(*outflowfsimap, false);

  // transposed fluid constraint matrix -> linearization of constraint equation
  *ConstrFluidMatrix = *FluidConstrMatrix->Transpose();
  ConstrFluidMatrix->Complete(fluidmap, constrmap);
  FluidConstrMatrix->Scale(invresscale);
  ConstrFluidMatrix->Scale(dttheta);

  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*DofRowMap(), true);
  outflowfsiinterface_->InsertCondVector(outflowfsiinterface_->ExtractCondVector(zeros), FluidRHS);
}


/*======================================================================*/
/* output of volume constraint related forces*/
void ADAPTER::FluidLung::OutputForces(Teuchos::RCP<Epetra_Vector> Forces)
{
  const Teuchos::RCP<IO::DiscretizationWriter>& output = DiscWriter();
  output->WriteVector("Add_Forces", Forces);
}
