/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.cpp

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "elch_algorithm.H"
// Output after each Outer Iteration step
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <iostream>
// fluid and transport solution are written to a file after each outer iteration loop
//#define OUTPUT

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false),
   natconv_(Teuchos::getIntegralValue<string>(prbdyn,"NATURAL_CONVECTION")),
   itmax_ (prbdyn.get<int>("ITEMAX")),
   ittol_ (prbdyn.get<double>("CONVTOL")),
   velincnp_ (rcp(new Epetra_Vector(*(FluidField().ExtractVelocityPart(FluidField().Velnp()))))),
   conpotincnp_ (rcp(new Epetra_Vector(*(ScaTraField().Phinp())))),
   density0_ (GetInitialFluidDensity())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoop()
{
  // write out inital state
  // Output();

  // provide information about initial state
  ScaTraField().OutputElectrodeInfo();
  ScaTraField().OutputMeanScalars();

  // ELCH algorithm without natural convection
  if (natconv_ == "No")
    TimeLoopElch();       // one-way coupling
  else
    TimeLoopConvection(); // two-way coupling (natural convection)

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoopElch()
{
    // compute error for problems with analytical solution
    ScaTraField().EvaluateErrorComparedToAnalyticalSol();

    // time loop
    while (NotFinished())
    {
      // prepare next time step
      PrepareTimeStep();

      // solve nonlinear Navier-Stokes system
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoopConvection()
{
    InitialCalculations();

    // time loop
    while(NotFinished())
    {
      IncrementTimeAndStep();

      // prepare next time step
      PrepareTimeStepConvection();

      // Outer iteration loop
      OuterIterationConvection();

      // update all single field solvers and density
      UpdateConvection();

      // write output to screen and files
      Output();
    }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::InitialCalculations()
{
  // a safety check
  switch((FluidField().TimIntScheme()))
  {
  case timeint_stationary:
  case timeint_one_step_theta:
  case timeint_bdf2:
    break;
  default: dserror("Selected time integration scheme is not available");
  }

  // compute initial density
  // set initial density to time step n+1 and n
  ScaTraField().ComputeDensity(density0_);
  // not essential
  ScaTraField().UpdateDensityElch();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  FluidField().PrepareTimeStep();

  // prepare time step
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::PrepareTimeStepConvection()
{
  /// LOMA functions are used for ELCH algorithm
  // only OST time integration scheme implemented:
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to
  // 0.0
  int numscal = 1;

  switch((FluidField().TimIntScheme()))
  {
  case timeint_stationary:
  case timeint_one_step_theta:
  case timeint_bdf2:
  {
    FluidField().SetIterLomaFields(
        ScaTraField().DensElchNp(),
        ScaTraField().DensElchN(),
        ScaTraField().DensElchNp(),
        1.0,
        1.0,
        0.0,
        numscal);
    break;
  }
  default: dserror("Selected time integration scheme is not available");
  }

  FluidField().PrepareTimeStep();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step()==1)
    ScaTraField().SetVelocityField(
        FluidField().Velnp(),
        FluidField().SgVelVisc(),
        Teuchos::null,
        FluidField().Discretization());

  // prepare time step (+ initialize one-step-theta scheme correctly with
  // velocity given above)
  ScaTraField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::DoFluidStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    FLUID SOLVER     \n";
    cout<<"*********************\n";
  }

  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  //FluidField().MultiCorrector();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::DoTransportStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    TRANSPORT SOLVER    \n";
    cout<<"************************\n";
  }

  // transfer actual velocity fields
  ScaTraField().SetVelocityField(
      FluidField().Velnp(),
      FluidField().SgVelVisc(),
      Teuchos::null,
      FluidField().Discretization()
  );

  // solve coupled transport equations for ion concentrations and
  // electric potential
  ScaTraField().NonlinearSolve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Update()
{
  FluidField().Update();

  // update scalar time derivative
  ScaTraField().UpdateTimeDerivative();

  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::UpdateConvection()
{
  ScaTraField().Update();

  // Old version:
  // ScaTraField().UpdateDensityElch();

  // LOMA functions are used for ELCH algorithm
  // only OST time integration scheme implemented:
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to
  // 0.0

  // Fluid field is solved once more to synchronizes the density and velocity field
  // Not a large influence since solution should be converged: Delta v in the range of 10E-6
  FluidField().NonlinearSolve();

  int numscal = 1;

  switch((FluidField().TimIntScheme()))
  {
  case timeint_stationary:
  case timeint_one_step_theta:
  case timeint_bdf2:
  {
    FluidField().SetIterLomaFields(
        ScaTraField().DensElchNp(),
        ScaTraField().DensElchN(),
        ScaTraField().DensElchNp(),
        1.0,
        1.0,
        0.0,
        numscal);
    break;
  }
  default: dserror("Selected time integration scheme is not available");
  }

  FluidField().Update();

  // Update density at time steps n+1 and n
  // Update density after SetTimeLomaFields
  ScaTraField().UpdateDensityElch();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::OuterIterationConvection()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  // Outer Iteration loop starts
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"**************************************************************\n";
    cout<<"      OUTER ONE-STEP-THETA/BDF2 ITERATION LOOP\n";
    printf("      Time Step %3d/%3d \n",Step(), ScaTraField().NStep());
    cout<<"**************************************************************\n";
  }

#ifdef OUTPUT
  // Output after each Outer Iteration step
  const int numdim = 3;
  //create output file name
  stringstream temp;
  temp << DRT::Problem::Instance()->OutputControlFile()->FileName()<<"_nonliniter_step"<<Step();
  string outname = temp.str();
  string probtype = DRT::Problem::Instance()->ProblemType(); // = "elch"

  RCP<IO::OutputControl> myoutputcontrol = rcp(new IO::OutputControl(ScaTraField().Discretization()->Comm(),probtype,"Polynomial","myinput",outname,numdim,0,1000));
  // create discretization writer with my own control settings
  RCP<IO::DiscretizationWriter> myoutput = rcp(new IO::DiscretizationWriter(ScaTraField().Discretization(),myoutputcontrol));
  // write mesh at step 0
  myoutput->WriteMesh(0,0.0);
#endif

  while (stopnonliniter==false)
  {
    itnum ++;

    conpotincnp_->Update(1.0,*ScaTraField().Phinp(),0.0);
    velincnp_->Update(1.0,*FluidField().ExtractVelocityPart(FluidField().Velnp()),0.0);

    // Density derivative is not used for OST, BDF2 and convective formulation
    int numscal = 1;
    if (itnum != 1)
    FluidField().SetTimeLomaFields(ScaTraField().DensElchNp(),0.0,null,numscal);

    // solve nonlinear Navier-Stokes system with body forces
    DoFluidStep();

    // solve nonlinear electrochemical transport equation
    DoTransportStep();
    //DoFluidStep();

    // compute new denselchnp_ and pass it to the fluid discretisation
    // LOMA functions are used in the ELCH algorithm
    // pass actual density field to fluid discretisation
    ScaTraField().ComputeDensity(density0_);

    // convergence check based on incremental values
    stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

    // Test: Output of the flux across the boundary into the output text file
    // after each outer Iteration step
    // print mean concentration
#ifdef OUTPUT
    if (stopnonliniter==false)
    {
    printf("\n");
    printf("Flux: Outer Iterations step: %3d \n", itnum);
    ScaTraField().OutputElectrodeInfo();
    ScaTraField().OutputMeanScalars();
    printf("\n");
    }

//#ifdef OUTPUT
   // iteration number (only after that data output is possible)
   myoutput->NewStep(itnum,itnum);
   myoutput->WriteVector("phinp", ScaTraField().Phinp());
   myoutput->WriteVector("convec_velocity", ScaTraField().ConVel());
   // myoutput->WriteVector("velnp", FluidField().Velnp());
#endif

  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool ELCH::Algorithm::ConvergenceCheck( int itnum,
    const int itmax,
    const double ittol)
{
  // convergence check based on the concentration, potential and
  // velocity increment

  //   | concentration increment |_2
  //  -------------------------------- < Tolerance
  //     | concentration_n+1 |_2

  bool stopnonliniter = false;
  LINALG::MapExtractor& conpotsplitter = ScaTraField().Splitter();
  // Variables to save different L2 - Norms

  double potincnorm_L2(0.0);
  double potnorm_L2(0.0);
  double velincnorm_L2(0.0);
  double velnorm_L2(0.0);
  double connorm_L2(0.0);
  double conincnorm_L2(0.0);

  // Calculate velocity increment and velocity L2 - Norm
  // velincnp_ = 1.0 * convelnp_ - 1.0 * conveln_

  velincnp_->Update(1.0,*FluidField().ExtractVelocityPart(FluidField().Velnp()),-1.0);
  velincnp_->Norm2(&velincnorm_L2); // Estimation of the L2 - norm save values to both variables (velincnorm_L2 and velnorm_L2)
  FluidField().ExtractVelocityPart(FluidField().Velnp())->Norm2(&velnorm_L2);

  // Calculate concentration increment and concentration L2 - Norm
  // tempincnp_ includes the concentration and the potential increment
  // tempincnp_ = 1.0 * phinp_ - 1.0 * phin_

  conpotincnp_->Update(1.0,*ScaTraField().Phinp(),-1.0);
  Teuchos::RCP<Epetra_Vector> onlycon = conpotsplitter.ExtractOtherVector(conpotincnp_);
  onlycon->Norm2(&conincnorm_L2);
  conpotsplitter.ExtractOtherVector(ScaTraField().Phinp(),onlycon);
  onlycon->Norm2(&connorm_L2);

  // Calculate potential increment and potential L2 - Norm

  Teuchos::RCP<Epetra_Vector> onlypot = conpotsplitter.ExtractCondVector(conpotincnp_);
  onlypot->Norm2(&potincnorm_L2);
  conpotsplitter.ExtractCondVector(ScaTraField().Phinp(),onlypot);
  onlypot->Norm2(&potnorm_L2);

  // care for the case that there is (almost) zero temperature or velocity
  // (usually not required for temperature)
  if (velnorm_L2 < 1e-6) velnorm_L2 = 1.0;
  if (connorm_L2 < 1e-6) connorm_L2 = 1.0;
  if (potnorm_L2 < 1e-6) potnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
    if (itnum != 1)
    {
      if (Comm().MyPID() == 0)
      {
      cout<<"\n";
      cout<<"*****************************************************************************\n";
      cout<<"                          OUTER ITERATION STEP\n";
      cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E  | %10.3E  |",
          itnum, itmax, ittol, conincnorm_L2/connorm_L2, potincnorm_L2/potnorm_L2, velincnorm_L2/velnorm_L2);
      printf("\n");
      printf("+------------+-------------------+-------------+--------------+-------------+\n");
      }

      // Converged or Not
      if ((conincnorm_L2/connorm_L2 <= ittol)  &&
          (potincnorm_L2/potnorm_L2 <= ittol) &&
          (velincnorm_L2/velnorm_L2 <= ittol))
        //if ((incconnorm_L2/connorm_L2 <= ittol))
      {
        stopnonliniter=true;
        if (Comm().MyPID() == 0)
        {
        printf("| Outer Iteration loop converged after iteration %3d/%3d                    |\n", itnum,itmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
        }
      }
      else
      {
        if (Comm().MyPID() == 0)
        {
        printf("| Outer Iteration loop is not converged after iteration %3d/%3d             |\n", itnum,itmax);
        printf("+---------------------------------------------------------------------- ----+");
        printf("\n");
        printf("\n");
        }
      }

      // warn if itemax is reached without convergence, but proceed to next timestep
      if (itnum == itmax)
      {
        if ((conincnorm_L2/connorm_L2 > ittol) ||
            (potincnorm_L2/potnorm_L2 > ittol) ||
            (velincnorm_L2/velnorm_L2 > ittol))
        {
          stopnonliniter=true;
          if ((Comm().MyPID() == 0))
          {
            printf("|     >>>>>> not converged in itemax steps!     |\n");
            printf("+-----------------------------------------------+\n");
            printf("\n");
            printf("\n");
          }
        }
      }
    }
    else
    {
      // first outer iteration loop: fluid solver has not got the new density yet
      // => minimum two outer iteration loops
      stopnonliniter=false;
      if (Comm().MyPID() == 0)
      {
      cout<<"\n";
      cout<<"*****************************************************************************\n";
      cout<<"                          OUTER ITERATION STEP\n";
      cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |       -      |      -      |      -      |",
          itnum, itmax, ittol);
      printf("\n");
      printf("+------------+-------------------+-------------+--------------+-------------+\n");
    }
  }

  return stopnonliniter;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ELCH::Algorithm::GetInitialFluidDensity()
{
  // initialization of the initial density
  double density = 0.0;

  // initialization only for convection, otherwise the density0 is set to zero
  if (natconv_ != "No")
  {
    // we ask the elements for the density (works for all fluid materials)
    ParameterList eleparams;
    eleparams.set("action","get_density");
    FluidField().Discretization()->Evaluate(eleparams,null,null,null,null,null);
    density = eleparams.get<double>("density");
    if (density <= 0.0) dserror("received illegal density value");
  }

  return density;
}

#endif // CCADISCRET
