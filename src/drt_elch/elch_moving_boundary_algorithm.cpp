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
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::MovingBoundaryAlgorithm::MovingBoundaryAlgorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidAleCouplingAlgorithm(comm,prbdyn,"FSICoupling"),
   molarvolume_(prbdyn.get<double>("MOLARVOLUME")),
   idispn_(FluidField().ExtractInterfaceVeln()),
   idispnp_(FluidField().ExtractInterfaceVeln()),
   iveln_(FluidField().ExtractInterfaceVeln())
{
  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);
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
  // write out initial state
  Output(); // do this to have ConditionID at all boundary conditions

  // compute error for problems with analytical solution
  ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    RCP<Epetra_Vector> incr(FluidField().ExtractInterfaceVeln());
    incr->PutScalar(0.0);
    double incnorm(0.0);
    int iter(0);

    // ToDo
    // improve this convergence test
    // (better check increment of ivel ????, test relative value etc.)
    while((iter == 0) or (incnorm > EPS10)) // do at least one step
    {
      iter++;

      /// compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_,iveln_);

      // save guessed value before solve
      incr->Update(1.0,*idispnp_,0.0);

      // solve nonlinear Navier-Stokes system on a deforming mesh
      SolveFluidAle();

      // solve transport equations for ion concentrations and electric potential
      SolveScaTra();

      /// compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_,iveln_);

      // compare with value after solving
      incr->Update(-1.0,*idispnp_,1.0);

      // compute L2 norm of increment
      incr->Norm2(&incnorm);

      if (Comm().MyPID()==0)
      {
        cout<<"After outer iteration "<<iter<<" ||dispincnp|| = "<<incnorm<<endl;
      }
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
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"  FLUID-ALE SOLVER   \n";
    cout<<"*********************\n";
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
void ELCH::MovingBoundaryAlgorithm::SolveFluidAle()
{
  //TEUCHOS_FUNC_TIME_MONITOR("ELCH::MovingBoundaryAlgorithm::SolveFluidAle");

  // solve nonlinear Navier-Stokes system on a moving mesh
  FluidAleNonlinearSolve(idispnp_,iveln_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::SolveScaTra()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"       ELCH SOLVER      \n";
    cout<<"************************\n";
  }

  if (FluidField().TimIntScheme()== INPAR::FLUID::timeint_gen_alpha)
  {
    dserror("ConvectiveVel() not implemented for Gen.Alpha");
  }
  else if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
  {
    dserror("ConvectiveVel() not implemented for AfGen.Alpha");
  }
  else
  {
    // transfer convective velocity = fluid velocity - grid velocity
    ScaTraField().SetVelocityField(
        FluidField().ConvectiveVel(), // = velnp - grid velocity
        FluidField().Hist(),
        Teuchos::null,
        FluidField().Discretization()
    );
  }

  // transfer moving mesh data
  ScaTraField().ApplyMeshMovement(
      FluidField().Dispnp(),
      FluidField().Discretization()
  );

  // solve coupled electrochemistry equations
  ScaTraField().Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Update()
{
  FluidField().Update();
  AleField().Update();
  ScaTraField().Update();

  // perform time shift of interface displacement
  idispn_->Update(1.0, *idispnp_ , 0.0);

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


void ELCH::MovingBoundaryAlgorithm::ComputeInterfaceVectors(
    RCP<Epetra_Vector> idispnp,
    RCP<Epetra_Vector> iveln)
{
  // calculate normal flux vector field only at FSICoupling boundaries (no output to file)
  vector<std::string> condnames(1);
  condnames[0] = "FSICoupling";
  Teuchos::RCP<Epetra_MultiVector> flux = ScaTraField().CalcFluxAtBoundary(condnames,false);

  // access discretizations
  RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  RCP<DRT::Discretization> scatradis = ScaTraField().Discretization();

  // no support for multiple reactions at the interface !
  // id of the reacting species
  int reactingspeciesid = 0;

  const Epetra_BlockMap& ivelmap = iveln->Map();

  // loop over all local nodes of fluid discretization
  for (int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); lnodeid++)
  {
    // Here we rely on the fact that the scatra discretization
    // is a clone of the fluid mesh. => a scatra node has the same
    // local (and global) ID as its corresponding fluid node!

    // get the processor's local fluid node with the same lnodeid
    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this fluid node
    vector<int> fluidnodedofs = fluiddis->Dof(fluidlnode);

    if(ivelmap.MyGID(fluidnodedofs[0])) // is this GID (implies: node) relevant for iveln_?
    {
      // determine number of space dimensions (numdof - pressure dof)
      const int numdim = ((int) fluidnodedofs.size()) -1;
      // number of dof per node in ScaTra
      int numscatradof = scatradis->NumDof(scatradis->lRowNode(lnodeid));

      vector<double> Values(numdim);
      for(int index=0;index<numdim;++index)
      {
        const int pos = lnodeid*numscatradof+reactingspeciesid;
        // interface growth has opposite direction of mass flow -> minus sign !!
        Values[index] = (-molarvolume_)*((*flux)[index])[pos];
      }

      // now insert only the first numdim entries (pressure dof is not inserted!)
      int error = iveln_->ReplaceGlobalValues(numdim,&Values[0],&fluidnodedofs[0]);
      if (error > 0) dserror("Could not insert values into vector iveln_: error %d",error);
    }
  }

  // have to compute an approximate displacement from given interface velocity
  // id^{n+1} = id^{n} + \delta t vel_i
  idispnp->Update(1.0,*idispn_,0.0);
  idispnp->Update(Dt(),*iveln_,1.0);

  return;
}


#endif // CCADISCRET
