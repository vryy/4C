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
   outmean_(Teuchos::getIntegralValue<int>(prbdyn,"OUTMEAN"))
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

    // solve nonlinear Navier-Stokes system on a deforming mesh
    DoFluidStep();

    // solve transport equations for ion concentrations and electric potential
    DoTransportStep();

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
//  Teuchos::RCP<Epetra_Vector> idisp_ = new FluidField().ExtractInterfaceVeln();
  const Teuchos::RCP<Epetra_Vector> idispn_ = rcp(new Epetra_Vector(*(FluidField().ExtractInterfaceVeln())));

  for (int lid = 0; lid < iveln_->MyLength(); ++lid)
  {
    if ((lid%3 == 0))
    {
      iveln_->ReplaceMyValue(lid+1,0,0.0025);
      double disp = 0.0025*Time();
      idispn_->ReplaceMyValue(lid+1,0,disp);
    }
    else
    {
      //idispn_->ReplaceMyValue(lid+1,0,0.0);
      //iveln_->ReplaceMyValue(lid+1,0,0.0);
    }
  }
  //cout << *idisp_;

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
      FluidField().Velnp(),
      FluidField().SgVelVisc(),
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
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().StatisticsAndOutput();
  ScaTraField().Output();
  if (outmean_) 
  {
    ScaTraField().OutputElectrodeInfo();
    ScaTraField().OutputMeanTempAndDens();
  }
  AleField().Output();

  return;
}


#endif // CCADISCRET
