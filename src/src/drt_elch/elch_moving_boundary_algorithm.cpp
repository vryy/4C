/*----------------------------------------------------------------------*/
/*!
\file elch_moving_boundary_algorithm.cpp

\brief Basis of all ELCH algorithms with moving boundaries

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "elch_moving_boundary_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::MovingBoundaryAlgorithm::MovingBoundaryAlgorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidAleCouplingAlgorithm(comm,prbdyn,"FSICoupling"),
   molarvolume_(prbdyn.get<double>("MOLARVOLUME")),
   idispn_(FluidField().ExtractInterfaceVeln())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::MovingBoundaryAlgorithm::~MovingBoundaryAlgorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::TimeLoop()
{
  // write out inital state
  // Output();

  // compute error for problems with analytical solution
  ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    for(int i = 0; i<2 ; i++)
    {
      if (Comm().MyPID()==0)
        cout<<"!!!!!!!!!! Sub-Iteration i = "<<i<<endl;

      // solve nonlinear Navier-Stokes system on a deforming mesh
      DoFluidStep();

      // solve transport equations for ion concentrations and electric potential
      DoTransportStep();
    }

    // update all single field solvers
    Update();

    // compute error for problems with analytical solution
    ScaTraField().EvaluateErrorComparedToAnalyticalSol();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n FLUID-ALE SOLVER \n******************\n";
  }

  FluidField().PrepareTimeStep();
  AleField().PrepareTimeStep();

  // prepare time step
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingMovingBoundaryAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::DoFluidStep()
{
  Teuchos::RCP<Epetra_Vector> iveln_ = FluidField().ExtractInterfaceVeln();

  // calculate normal flux vector field only at FSICoupling boundaries
  vector<std::string> condnames;
  condnames.push_back("FSICoupling");
  condnames.push_back("ElectrodeKinetics");
  Teuchos::RCP<Epetra_MultiVector> flux = ScaTraField().CalcFluxAtBoundary(condnames);

  RCP<DRT::Discretization> fluiddis = FluidField().Discretization();

  // number of dof per node in ScaTra
  RCP<DRT::Discretization> scatradis = ScaTraField().Discretization();
  int numscatradof = scatradis->NumDof(scatradis->lRowNode(0));

  // id of the reacting species
  int reactingspeciesid = 0; // default

  // loop over all local nodes of scatra discretization
  for (int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); lnodeid++)
  {
    // Here we rely on the fact that the scatra discretization
    // is a clone of the fluid mesh. => a scatra node has the same
    // local (and global) ID as its corresponding fluid node!

    // get the processor's local fluid node with the same lnodeid
    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this fluid node
    vector<int> fluidnodedofs = fluiddis->Dof(fluidlnode);
    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;

    vector<double> Values(numdim);

    for(int index=0;index<numdim;++index)
    {
      const int pos = lnodeid*numscatradof+reactingspeciesid;
      // interface growth has opposite direction of mass flow -> minus sign
      Values[index] = (-molarvolume_)*((*flux)[index])[pos];
    }

    // now insert only the first numdim entries (nothing will happen for
    // nodes not belonging to the FSI-Coupling -> error == 1)
 //int error = iveln_->ReplaceGlobalValues(numdim,&Values[0],&fluidnodedofs[0]);
    iveln_->ReplaceGlobalValues(numdim,&Values[0],&fluidnodedofs[0]);
  }

  // have to compute the approximate displacement from interface velocity
  double timescale = 1./FluidField().TimeScaling();
  idispn_->Update(Dt(),*iveln_,timescale);

  iveln_->PutScalar(0.0);

  // solve nonlinear Navier-Stokes system on a moving mesh
  FluidAleNonlinearSolve(idispn_,iveln_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::DoTransportStep()
{

  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n TRANSPORT SOLVER \n******************\n";
  }

  // transfer actual velocity fields
  ScaTraField().SetVelocityField(
      FluidField().ConvectiveVel(), // = velnp - grid velocity
      FluidField().Hist(),
      Teuchos::null,
      FluidField().Discretization()
  );

  // transfer moving mesh data
  ScaTraField().ApplyMeshMovement(
      FluidField().Dispnp(),
      FluidField().Discretization()
  );

  // solve coupled transport equations for ion concentrations and electric
  // potential
  ScaTraField().NonlinearSolve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Update()
{
  FluidField().Update();
  AleField().Update();
  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();
  ScaTraField().Output();
  AleField().Output();

  return;
}


#endif // CCADISCRET
