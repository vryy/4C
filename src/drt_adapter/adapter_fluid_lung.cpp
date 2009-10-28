/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
#include "adapter_fluid_lung.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"

/*======================================================================*/
/* constructor */
ADAPTER::FluidLung::FluidLung(Teuchos::RCP<Fluid> fluid)
: FluidWrapper(fluid)
{
  // make sure
  if (fluid_ == null)
    dserror("Failed to create the underlying fluid adapter");

  // get lung fluid-structure volume constraints
  std::vector<DRT::Condition*> temp;
  Discretization()->GetCondition("StructFluidSurfCoupling",temp);
  for (unsigned i=0; i<temp.size(); ++i)
  {
    DRT::Condition& cond = *(temp[i]);
    if (*(cond.Get<string>("field")) == "fluid")
      constrcond_.push_back(temp[i]);
  }
  if (constrcond_.size() == 0) dserror("No structure-fluid volume constraints found for lung fsi");
}


/*======================================================================*/
/* list of coupled fluid-structure volumes */
void ADAPTER::FluidLung::ListLungVolCons(std::set<int>& LungVolConIDs,
                                         int& MinLungVolConID)
{
  MinLungVolConID = 1;

  for (unsigned int i=0; i<constrcond_.size();++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);
    int condID = cond.GetInt("coupling id");
    if (condID < MinLungVolConID)
      MinLungVolConID = condID;
    LungVolConIDs.insert(condID);
  }
}


/*======================================================================*/
/* determine initial flow rates */
void ADAPTER::FluidLung::InitializeVolCon(Teuchos::RCP<Epetra_Vector> initflowrate,
                                          const int offsetID)
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

    //if (*(cond.Get<string>("field")) == "fluid")
    {
      // Get ConditionID of current condition if defined and write value in parameterlist

      int condID=cond.GetInt("coupling id");

      ParameterList params;
      params.set("ConditionID",condID);
      params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));
      params.set("action","flowrate_deriv");
      params.set("flowrateonly", true);
      const double dt = Dt();
      params.set("dt", dt);

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*Discretization(),lm,lmowner);

        // Reshape element matrices and vectors and init to zero
        elevector3.Size(1);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*Discretization(),lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly

        vector<int> constrlm;
        vector<int> constrowner;
        constrlm.push_back(condID-offsetID);
        constrowner.push_back(curr->second->Owner());
        LINALG::Assemble(*initflowrate,elevector3,constrlm,constrowner);
      }
    }
  }
}


/*======================================================================*/
/* evaluate structural part of fluid-structure volume constraint */
void ADAPTER::FluidLung::EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> FluidShapeDerivMatrix,
                                        Teuchos::RCP<LINALG::SparseMatrix> FluidConstrMatrix,
                                        Teuchos::RCP<LINALG::SparseMatrix> ConstrFluidMatrix,
                                        Teuchos::RCP<LINALG::BlockSparseMatrixBase> AleConstrMatrix,
                                        Teuchos::RCP<LINALG::BlockSparseMatrixBase> ConstrAleMatrix,
                                        Teuchos::RCP<Epetra_Vector> FluidRHS,
                                        Teuchos::RCP<Epetra_Vector> CurrFlowRates,
                                        Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
                                        const int offsetID,
                                        const double dttheta)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // set ale displacements, fluid and grid velocities
  Discretization()->ClearState();
  Discretization()->SetState("convectivevel", ConvectiveVel());
  Discretization()->SetState("dispnp", Dispnp());

  // residual scaling for fluid matrices needed
  const double invresscale = 1.0/ResidualScaling();

  //---------------------------------------------------------------------
  // loop through conditions and evaluate them
  //---------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    //if (*(cond.Get<string>("field")) == "fluid")
    {
      // Get ConditionID of current condition if defined and write value in parameterlist
      int condID=cond.GetInt("coupling id");
      ParameterList params;
      params.set("ConditionID",condID);
      const double dt = Dt();
      params.set("dt", dt);

      // global and local ID of this bc in the redundant vectors
      int gindex = condID-offsetID;
      const int lindex = (CurrFlowRates->Map()).LID(gindex);

      // Get the current lagrange multiplier value for this condition
      const double lagraval = (*lagrMultVecRed)[lindex];

      // elements might need condition
      params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;  // (d^2 Q)/(du dd)
      Epetra_SerialDenseMatrix elematrix2;  // (d^2 Q)/(dd)^2
      Epetra_SerialDenseVector elevector1;  // dQ/du
      Epetra_SerialDenseVector elevector2;  // dQ/dd
      Epetra_SerialDenseVector elevector3;  // Q

      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // no check for empty geometry here since in parallel computations
      // there might be processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;

      // define element action
      params.set("action", "flowrate_deriv");

      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*Discretization(),lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.Shape(eledim,eledim);
        elematrix2.Shape(eledim,eledim);
        elevector1.Size(eledim);
        elevector2.Size(eledim);
        elevector3.Size(1);

        //---------------------------------------------------------------------
        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*Discretization(),lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        //---------------------------------------------------------------------
        // assembly
        int eid = curr->second->Id();

        elematrix1.Scale(-lagraval*invresscale*dttheta);
        FluidShapeDerivMatrix->Assemble(eid,elematrix1,lm,lmowner);

        // assemble to rectangular matrix. The column corresponds to the constraint ID.
        vector<int> colvec(1);
        colvec[0]=gindex;

        elevector1.Scale(-dttheta*invresscale);
        FluidConstrMatrix->Assemble(eid,elevector1,lm,lmowner,colvec);

        elevector2.Scale(-dttheta);
        AleConstrMatrix->Assemble(eid,elevector2,lm,lmowner,colvec);

        // negative sign (for shift to rhs) is already implicitly taken into account!
        elevector1.Scale(-lagraval);
        LINALG::Assemble(*FluidRHS,elevector1,lm,lmowner);

        vector<int> constrlm;
        vector<int> constrowner;
        constrlm.push_back(gindex);
        constrowner.push_back(curr->second->Owner());
        LINALG::Assemble(*CurrFlowRates,elevector3,constrlm,constrowner);
      }
    }
  }

  FluidShapeDerivMatrix->Complete();

  const Epetra_Map& constrmap = ConstrFluidMatrix->RangeMap();
  const Epetra_Map& fluidmap = FluidConstrMatrix->RangeMap();
  FluidConstrMatrix->Complete(constrmap, fluidmap);
  AleConstrMatrix->Complete();

  // note: since FluidConstrMatrix was multiplied by invresscale, we
  // have to divide by it now for obtaining ConstrFluidMatrix
  ConstrFluidMatrix->Add(*FluidConstrMatrix, true, 1.0/invresscale, 0.0);
  ConstrFluidMatrix->Complete(fluidmap, constrmap);

  for (int i=0; i<Interface().NumMaps(); ++i)
    ConstrAleMatrix->Matrix(0,i).Add(AleConstrMatrix->Matrix(i,0), true, 1.0, 0.0);
  ConstrAleMatrix->Complete();

  // Apply Dirichlet BC to relevant stiffness matrices
  FluidShapeDerivMatrix->ApplyDirichlet(*GetDBCMapExtractor()->CondMap(), false);
  FluidConstrMatrix->ApplyDirichlet(*GetDBCMapExtractor()->CondMap(), false);
  AleConstrMatrix->ApplyDirichlet(*GetDBCMapExtractor()->CondMap(), false);

  // Blank residual at DOFs on Dirichlet BC
  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*DofRowMap(), true);
  GetDBCMapExtractor()->InsertCondVector(GetDBCMapExtractor()->ExtractCondVector(zeros), FluidRHS);
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
