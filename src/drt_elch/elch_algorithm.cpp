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

#include "elch_algorithm.H"
#include "../drt_scatra/scatra_utils.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../linalg/linalg_mapextractor.H"
// Output after each Outer Iteration step
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <iostream>
// fluid and transport solution are written to a file after each outer iteration loop
//#define OUTPUT

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,"scatra",solverparams),
   natconv_(DRT::INPUT::IntegralValue<int>(prbdyn,"NATURAL_CONVECTION")),
   itmax_ (prbdyn.get<int>("ITEMAX")),
   ittol_ (prbdyn.get<double>("CONVTOL")),
   velincnp_ (Teuchos::rcp(new Epetra_Vector(*(FluidField().ExtractVelocityPart(FluidField().Velnp()))))),
   conpotincnp_ (Teuchos::rcp(new Epetra_Vector(*(ScaTraField().Phinp())))),
   samstart_(prbdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_START")),
   samstop_(prbdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_STOP"))
{
  // Setup of TurbulenceStatisticManager is performed in the
  // constructor of class ScaTraFluidCouplingAlgorithm

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
  // provide information about initial field (do not do for restarts!)
  if (Step()==0)
  {
    // write out initial state
    // Output();

    ScaTraField().OutputElectrodeInfo();
    ScaTraField().OutputMeanScalars();

    // compute error for problems with analytical solution (initial field!)
    ScaTraField().EvaluateErrorComparedToAnalyticalSol();
  }

  // switch ELCH algorithm
  if (natconv_ == false)
    TimeLoopElch();       // one-way coupling
  else
    TimeLoopConvection(); // two-way coupling (natural convection)

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoopElch()
{
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
  case INPAR::FLUID::timeint_stationary:
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_bdf2:
    break;
  default: dserror("Selected time integration scheme is not available");
    break;
  }

  // compute initial density
  // set initial density to time step n+1 and n
  ScaTraField().ComputeDensity();
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
  // OST/BDF2 time integration schemes are implemented
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST and BDF2 (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to 0.0

  switch((FluidField().TimIntScheme()))
  {
  case INPAR::FLUID::timeint_stationary:
  {
    FluidField().SetIterLomaFields(
        ScaTraField().DensElchNp(),
        ScaTraField().DensElchNp(), // we have to provide something here
        Teuchos::null,
        Teuchos::null,
        1.0,
        1.0,
        0.0,
        0.0,
        ScaTraField().Discretization());
    break;
  }
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_bdf2:
  {
    FluidField().SetIterLomaFields(
        ScaTraField().DensElchNp(),
        ScaTraField().DensElchN(),
        Teuchos::null,
        Teuchos::null,
        1.0,
        1.0,
        0.0,
        0.0,
        ScaTraField().Discretization());
    break;
  }
  default: dserror("Selected time integration scheme is not available");
    break;
  }

  FluidField().PrepareTimeStep();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step()==1)
    ScaTraField().SetVelocityField(
        FluidField().Velnp(),
        FluidField().Hist(),
        Teuchos::null,
        Teuchos::null,
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
    std::cout<<"\n";
    std::cout<<"*********************\n";
    std::cout<<"    FLUID SOLVER     \n";
    std::cout<<"*********************\n";
  }

  // solve nonlinear Navier-Stokes system
  FluidField().Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::DoTransportStep()
{
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n";
    std::cout<<"************************\n";
    std::cout<<"       ELCH SOLVER      \n";
    std::cout<<"************************\n";
  }

  // transfer convective velocity to scalar transport field solver
  switch(FluidField().TimIntScheme())
  {
  case INPAR::FLUID::timeint_npgenalpha:
  case INPAR::FLUID::timeint_afgenalpha:
  {
    ScaTraField().SetVelocityField(
        FluidField().Velaf(),
        FluidField().Accam(),
        Teuchos::null,
        Teuchos::null,
        Teuchos::null,
        FluidField().Discretization());
  }
  break;
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_bdf2:
  case INPAR::FLUID::timeint_stationary:
  {
    ScaTraField().SetVelocityField(
        FluidField().Velnp(),
        FluidField().Hist(),
        Teuchos::null,
        Teuchos::null,
        Teuchos::null,
        FluidField().Discretization()
    );
  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }

  // solve coupled transport equations for ion concentrations and
  // electric potential
  ScaTraField().Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::UpdateConvection()
{
  // Update ScaTra fields
  ScaTraField().Update();

  // OST/BDF2 time integration schemes are implemented
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST and BDF2 (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to 0.0

  // here SetIterLomaFields is only necessary if the density is used to update the accelerations
  // otherwise it is a redundant step since the density is multiplied to the history vector in the element
  /*
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
        Teuchos::null,
        Teuchos::null,
        1.0,
        1.0,
        0.0,
        0.0,
        numscal);
    break;
  }
  default: dserror("Selected time integration scheme is not available");
  }
  */

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

  if ((Step()>=samstart_) and (Step()<=samstop_))
  {
    // if statistics for one-way coupled problems is performed, provide
    // the field for the first scalar!
    FluidField().SetTimeLomaFields(
        ScaTraField().Phinp(),
        0.0,
        Teuchos::null,
        ScaTraField().Discretization(),
        0 // do statistics for FIRST dof at every node!!
    );
  }

  FluidField().StatisticsAndOutput();
  ScaTraField().Output();

  // we have to call the output of averaged fields for scatra separately
  if (  FluidField().TurbulenceStatisticManager() != Teuchos::null)
    FluidField().TurbulenceStatisticManager()
        ->DoOutputForScaTra(ScaTraField().DiscWriter(),ScaTraField().Step());

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
    std::cout<<"\n";
    std::cout<<"**************************************************************\n";
    std::cout<<"      OUTER ITERATION LOOP ("<<ScaTraField().MethodTitle()<<")\n";
    printf("      Time Step %3d/%3d \n",Step(), ScaTraField().NStep());
    std::cout<<"**************************************************************\n";
  }

#ifdef OUTPUT
  // Output after each Outer Iteration step
  const int numdim = 3;
  //create output file name
  std::stringstream temp;
  temp << DRT::Problem::Instance()->OutputControlFile()->FileName()<<"_nonliniter_step"<<Step();
  std::string outname = temp.str();
  std::string probtype = DRT::Problem::Instance()->ProblemName();

  RCP<IO::OutputControl> myoutputcontrol = Teuchos::rcp(new IO::OutputControl(ScaTraField().Discretization()->Comm(),probtype,"Polynomial","myinput",outname,numdim,0,1000));
  // create discretization writer with my own control settings
  RCP<IO::DiscretizationWriter> myoutput = Teuchos::rcp(new IO::DiscretizationWriter(ScaTraField().Discretization(),myoutputcontrol));
  // write mesh at step 0
  myoutput->WriteMesh(0,0.0);
#endif

  while (stopnonliniter==false)
  {
    itnum ++;

    conpotincnp_->Update(1.0,*ScaTraField().Phinp(),0.0);
    velincnp_->Update(1.0,*FluidField().ExtractVelocityPart(FluidField().Velnp()),0.0);

    // solve nonlinear Navier-Stokes system with body forces
    DoFluidStep();

    // solve nonlinear electrochemical transport equation
    DoTransportStep();

    // compute new denselchnp_ and pass it to the fluid discretisation
    // pass actual density field to fluid discretisation
    // Density derivative is not used for OST, BDF2 and convective formulation
    ScaTraField().ComputeDensity();
    FluidField().SetTimeLomaFields(
        ScaTraField().DensElchNp(),
        0.0,
        Teuchos::null,
        ScaTraField().Discretization());

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
  Teuchos::RCP<LINALG::MapExtractor> conpotsplitter = ScaTraField().Splitter();
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
  if (SCATRA::IsElch(ScaTraField().ScaTraType()))
  {
  Teuchos::RCP<Epetra_Vector> onlycon = conpotsplitter->ExtractOtherVector(conpotincnp_);
  onlycon->Norm2(&conincnorm_L2);
  conpotsplitter->ExtractOtherVector(ScaTraField().Phinp(),onlycon);
  onlycon->Norm2(&connorm_L2);

  // Calculate potential increment and potential L2 - Norm
  Teuchos::RCP<Epetra_Vector> onlypot = conpotsplitter->ExtractCondVector(conpotincnp_);
  onlypot->Norm2(&potincnorm_L2);
  conpotsplitter->ExtractCondVector(ScaTraField().Phinp(),onlypot);
  onlypot->Norm2(&potnorm_L2);
  }
  else // only one scalar present (convection-diffusion)
  {
    conpotincnp_->Norm2(&conincnorm_L2);
    ScaTraField().Phinp()->Norm2(&connorm_L2);
  }

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
      std::cout<<"\n";
      std::cout<<"*****************************************************************************\n";
      std::cout<<"                          OUTER ITERATION STEP\n";
      std::cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E  | %10.3E  |",
          itnum, itmax, ittol, conincnorm_L2/connorm_L2, potincnorm_L2/potnorm_L2, velincnorm_L2/velnorm_L2);
      printf("\n");
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
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
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
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
      std::cout<<"\n";
      std::cout<<"*****************************************************************************\n";
      std::cout<<"                          OUTER ITERATION STEP\n";
      std::cout<<"*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |       -      |      -      |      -      |",
          itnum, itmax, ittol);
      printf("\n");
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
    }
  }

    // warn if itemax is reached without convergence, but proceed to next timestep
    // itemax = 1 is also possible for segregated coupling approaches (not fully implicit)
    if (itnum == itmax)
    {
      if (((conincnorm_L2/connorm_L2 > ittol) ||
          (potincnorm_L2/potnorm_L2 > ittol) ||
          (velincnorm_L2/velnorm_L2 > ittol)) or (itmax==1))
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

  return stopnonliniter;
}
