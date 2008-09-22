/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.cpp

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "loma_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::Algorithm(
    Epetra_Comm& comm,
    RCP<DRT::Discretization>      scatradis,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn),
   scatradiscret_(scatradis)
{
  // get noderowmap of scatra discretization
  const Epetra_Map* scatranodrowmap = scatradiscret_->NodeRowMap();

  // node-based densities at time steps n+1, n and n-1
  noddensnp_ = LINALG::CreateVector(*scatranodrowmap,true);
  noddensn_  = LINALG::CreateVector(*scatranodrowmap,true);
  noddensnm_ = LINALG::CreateVector(*scatranodrowmap,true);

  // node-based density increment for convergence check
  densincnp_ = LINALG::CreateVector(*scatranodrowmap,true);

  // node-based temperatures at time step n+1
  nodtempnp_ = LINALG::CreateVector(*scatranodrowmap,true);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeLoop()
{
  // initial thermodynamic pressure (in N/m� = kg/(m*s�) = J/m�)
  // constantly set to atmospheric pressure, for the time being -> dp_therm/dt=0
  thermpressnp_ = 98100.0;

  // get initial temperature field for computation of initial density field
  GetCurrentScalar();

  // compute initial density field using initial temperature + therm. pressure
  ComputeDensity();
  noddensn_->Update(1.0,*noddensnp_,0.0);

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // predict density field (from second time step on)
    if (NotFirstTimeStep())
    {
      // predictor: linear extrapolation in time
      //noddensnp_->Update(2.0,*noddensn_,-1.0,*noddensnm_,0.0);

      // predictor: constant extrapolation in time
      noddensnp_->Update(1.0,*noddensn_,0.0);
    }

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop
    OuterLoop();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // set field vectors: initial density-weighted convective velocity + density
  //if (Step()==1) ScaTraField().SetIterLomaFields(2,ConvectiveVelocity(),noddensn_);

  // set field vectors: density at n and n-1
  ScaTraField().SetTimeLomaFields(noddensn_,noddensnm_);

  // prepare temperature time step (+ initialize one-step-theta scheme correctly)
  ScaTraField().PrepareTimeStep();

  // set field vectors: density at n and n-1
  FluidField().SetTimeLomaFields(noddensn_,noddensnm_);

  // prepare fluid time step, particularly predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OuterLoop()
{
  const double  ittol = 1e-3;
  int  itnum = 0;
  int  itemax = 2;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
    cout<<"\n**********************\n OUTER ITERATION LOOP \n**********************\n";

  while (stopnonliniter==false)
  {
    itnum++;

    // get new velocity field
    GetCurrentFluidVelocity();

    // set field vectors: density-weighted convective velocity + density
    ScaTraField().SetIterLomaFields(2,ConvectiveVelocity(),noddensnp_);

    // solve transport equation for temperature
    if (Comm().MyPID()==0)
      cout<<"\n**********************\n  TEMPERATURE SOLVER \n**********************\n";
    ScaTraField().Solve();

    // get temperature
    GetCurrentScalar();

    // compute density using current temperature + thermodynamic pressure
    ComputeDensity();

    // set field vectors: density
    FluidField().SetIterLomaFields(noddensnp_);

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0)
      cout<<"\n*********************\n     FLOW SOLVER \n*********************\n";
    FluidField().NonlinearSolve();

    double densincnorm_L2;
    double densnorm_L2;

    densincnp_->Update(1.0,*noddensnp_,-1.0,*noddensn_,0.0);
    densincnp_->Norm2(&densincnorm_L2);
    noddensnp_->Norm2(&densnorm_L2);

    // care for the case that there is (almost) zero density
    if (densnorm_L2 < 1e-6) densnorm_L2 = 1.0;

    if (Comm().MyPID()==0)
    {
      cout<<"\n**********************\n OUTER ITERATION STEP \n**********************\n";
      printf("+------------+-------------------+--------------+--------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- dens-norm -|-- dens-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
             itnum,itemax,ittol,densnorm_L2,densincnorm_L2/densnorm_L2);
      printf("\n");
      printf("+------------+-------------------+--------------+--------------+\n");
    }

    if (densincnorm_L2/densnorm_L2 <= ittol) stopnonliniter=true;

    // warn if itemax is reached without convergence, but proceed to next timestep
    if ((itnum == itemax) and (densincnorm_L2/densnorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (Comm().MyPID()==0)
      {
        printf("|    >>>>>> not converged in itemax steps!                     |\n");
        printf("+--------------------------------------------------------------+\n");
      }
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::ComputeDensity()
{
  // set specific gas constant (in J/(kg*K))
  const double gasconst = 287.05;

  // set node-based temperature field vector out of dof-based vector
  nodtempnp_->Update(1.0,*Scalar(),0.0);

  // compute density based on equation of state
  double fac = thermpressnp_/gasconst;
  noddensnp_->Reciprocal(*nodtempnp_);
  noddensnp_->Scale(fac);

  //noddensnp_->PutScalar(1.0);
  //noddensn_->PutScalar(1.0);
  //noddensnm_->PutScalar(1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();

  // update of node-based density solutions
  noddensnm_->Update(1.0,*noddensn_ ,0.0);
  noddensn_->Update(1.0,*noddensnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Output()
{
  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
