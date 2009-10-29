/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
#include "adapter_structure_lung.H"
#include "../drt_io/io.H"

/*======================================================================*/
/* constructor */
ADAPTER::StructureLung::StructureLung(Teuchos::RCP<Structure> stru)
: StructureWrapper(stru)
{
  // make sure
  if (structure_ == null)
    dserror("Failed to create the underlying structural adapter");

  // get lung fluid-structure volume constraints
  std::vector<DRT::Condition*> temp;
  Discretization()->GetCondition("StructFluidSurfCoupling",temp);
  for (unsigned i=0; i<temp.size(); ++i)
  {
    DRT::Condition& cond = *(temp[i]);
    if (*(cond.Get<string>("field")) == "structure")
      constrcond_.push_back(temp[i]);
  }
  if (constrcond_.size() == 0) dserror("No structure-fluid volume constraints found for lung fsi");

  // build mapextractor for fsi <-> full map
  fsiinterface_ = LINALG::MapExtractor(*Interface().FullMap(), Interface().FSICondMap());
}


/*======================================================================*/
/* list of coupled fluid-structure volumes */
void ADAPTER::StructureLung::ListLungVolCons(std::set<int>& LungVolConIDs,
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
/* determine initial volumes */
void ADAPTER::StructureLung::InitializeVolCon(Teuchos::RCP<Epetra_Vector> initvol,
                                              Teuchos::RCP<Epetra_Vector> signvol,
                                              const int offsetID)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  if (initvol == Teuchos::null or signvol == Teuchos::null)
    dserror("Structure lung volume constraint cannot be initialized");

  if (!initvol->Map().SameAs(signvol->Map()))
    dserror("Structure lung volume constraint cannot be initialized");

  ParameterList params;
  params.set("action","calc_struct_constrvol");

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
      int condID=cond.GetInt("coupling id");
      params.set("ConditionID",condID);
      params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));

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

        // reshape element matrices and vectors and init to zero
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
        LINALG::Assemble(*initvol,elevector3,constrlm,constrowner);
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

  for (int i=0; i<initvol->MyLength(); ++i)
  {
    if ((*initvol)[i] > 0.0)
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
void ADAPTER::StructureLung::EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> StructMatrix,
                                            Teuchos::RCP<Epetra_Vector> StructRHS,
                                            Teuchos::RCP<Epetra_Vector> CurrVols,
                                            Teuchos::RCP<Epetra_Vector> SignVols,
                                            Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
                                            const int offsetID)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // Define element action
  ParameterList params;
  params.set("action", "calc_struct_volconstrstiff");

  // set displacements
  Discretization()->ClearState();
  Discretization()->SetState("displacement", Dispnp());

  //---------------------------------------------------------------------
  // loop through conditions and evaluate them
  //---------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    //cout << "CONDITION " << i << endl;

    DRT::Condition& cond = *(constrcond_[i]);

    //if (*(cond.Get<string>("field")) == "structure")
    {
      // Get ConditionID of current condition if defined and write value in parameterlist
      int condID = cond.GetInt("coupling id");
      params.set("ConditionID",condID);

      // global and local ID of this bc in the redundant vectors
      int gindex = condID-offsetID;
      const int lindex = (CurrVols->Map()).LID(gindex);
      if (lindex == -1)
        dserror("Corrupt vector of current volumes");

      // Get the current lagrange multiplier value for this condition
      const double lagraval = (*lagrMultVecRed)[lindex];

      // Get the sign (of volume) for this condition
      const double sign = (*SignVols)[lindex];

      // elements might need condition
      params.set<RefCountPtr<DRT::Condition> >("condition", rcp(&cond,false));

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // no check for empty geometry here since in parallel computations
      // there might be processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;

      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*Discretization(),lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.Shape(eledim,eledim);  // stiffness part
        elematrix2.Shape(eledim,eledim);  // this element matrix is only needed for the function call only
        elevector1.Size(eledim);          // rhs part
        elevector2.Size(eledim);          // constraint matrix
        elevector3.Size(1);               // current volume

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*Discretization(),lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly
        int eid = curr->second->Id();

        // NOTE:
        // - no time integration related scaling needed here, since everything is evaluated at the
        //   end of the time step and assembled to the already completed system.
        // - all element quantities are multiplied by the previously derived sign in order to assure
        //   that computed volumes are positive and corresponding derivatives are determined with a
        //   consistent sign.
        // - additional multiplication with -1.0 accounts for different definition of volume
        //   difference ("normal" volume constraint: Vref-Vcurr, lung volume constraint: Vcurr-Vref,
        //   which seems more natural in this case)

        elematrix1.Scale(-lagraval*sign);
        StructMatrix->Assemble(eid,elematrix1,lm,lmowner);

        // assemble to rectangular matrix. The column corresponds to the constraint ID.
        vector<int> colvec(1);
        colvec[0]=gindex;
        elevector2.Scale(-sign);
        StructMatrix->Assemble(eid,elevector2,lm,lmowner,colvec);

        // "Newton-ready" residual -> already scaled with -1.0
        elevector1.Scale(lagraval*sign);
        LINALG::Assemble(*StructRHS,elevector1,lm,lmowner);

        // No scaling with -1.0 necessary here, since the constraint rhs is determined consistently,
        // i.e.  -(Vcurr - Vold) in the fsi algorithm, thus -1.0 is included there.
        vector<int> constrlm;
        vector<int> constrowner;
        constrlm.push_back(gindex);
        constrowner.push_back(curr->second->Owner());
        elevector3.Scale(sign);
        LINALG::Assemble(*CurrVols,elevector3,constrlm,constrowner);
      }
    }
  }

  StructMatrix->Complete();

  const Teuchos::RCP<const Epetra_Map >& dispmap = StructMatrix->RangeExtractor().Map(0);

  LINALG::SparseMatrix& ConstrStructMatrix = StructMatrix->Matrix(1,0);
  ConstrStructMatrix.UnComplete();
  ConstrStructMatrix.Add(StructMatrix->Matrix(0,1), true, 1.0, 0.0);

  StructMatrix->Complete();

  // Apply Dirichlet BC to stiffness matrix

  StructMatrix->ApplyDirichlet(*GetDBCMapExtractor()->CondMap(), false);

  // blank residual at DOFs on Dirichlet BC
  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dispmap, true);
  GetDBCMapExtractor()->InsertCondVector(GetDBCMapExtractor()->ExtractCondVector(zeros), StructRHS);

}


/*======================================================================*/
/* write additional stuff in case of restart */
void ADAPTER::StructureLung::WriteVolConRestart(Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
                                                Teuchos::RCP<Epetra_Vector> OldVolsRed,
                                                Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
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

  Teuchos::RCP<std::vector<double> > flowrates = rcp(new std::vector<double>(OldFlowRatesRed->MyLength()));
  for (int i=0; i<OldFlowRatesRed->MyLength(); ++i)
  {
    (*flowrates)[i] = (*OldFlowRatesRed)[i];
  }
  Teuchos::RCP<std::vector<double> > volumes = rcp(new std::vector<double>(OldVolsRed->MyLength()));
  for (int i=0; i<OldVolsRed->MyLength(); ++i)
  {
    (*volumes)[i] = (*OldVolsRed)[i];
  }
  Teuchos::RCP<std::vector<double> > lmult = rcp(new std::vector<double>(OldLagrMultRed->MyLength()));
  for (int i=0; i<OldLagrMultRed->MyLength(); ++i)
  {
    (*lmult)[i] = (*OldLagrMultRed)[i];
  }

  output->WriteRedundantDoubleVector(stream1.str(), volumes);
  output->WriteRedundantDoubleVector(stream2.str(), flowrates);
  output->WriteRedundantDoubleVector(stream3.str(), lmult);
}


/*======================================================================*/
/* read additional stuff in case of restart */
void ADAPTER::StructureLung::ReadVolConRestart(const int step,
                                               Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
                                               Teuchos::RCP<Epetra_Vector> OldVolsRed,
                                               Teuchos::RCP<Epetra_Vector> OldLagrMultRed)
{
  IO::DiscretizationReader reader(Discretization(), step);
  if (step != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");
  std::stringstream stream1;
  stream1 << "OldVols";
  std::stringstream stream2;
  stream2 << "OldFlowRates";
  std::stringstream stream3;
  stream3 << "OldLagrMult";

  Teuchos::RCP<std::vector<double> > flowrates = rcp(new std::vector<double>(OldFlowRatesRed->MyLength()));
  Teuchos::RCP<std::vector<double> > volumes = rcp(new std::vector<double>(OldVolsRed->MyLength()));
  Teuchos::RCP<std::vector<double> > lmult = rcp(new std::vector<double>(OldLagrMultRed->MyLength()));

  reader.ReadRedundantDoubleVector(volumes, stream1.str());
  reader.ReadRedundantDoubleVector(flowrates, stream2.str());
  reader.ReadRedundantDoubleVector(lmult, stream3.str());

  for (int i=0; i<OldFlowRatesRed->MyLength(); ++i)
  {
    (*OldFlowRatesRed)[i] = (*flowrates)[i];
  }
  for (int i=0; i<OldVolsRed->MyLength(); ++i)
  {
    (*OldVolsRed)[i] = (*volumes)[i];
  }
  for (int i=0; i<OldLagrMultRed->MyLength(); ++i)
  {
    (*OldLagrMultRed)[i] = (*lmult)[i];
  }
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
