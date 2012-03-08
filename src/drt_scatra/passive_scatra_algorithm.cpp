/*----------------------------------------------------------------------*/
/*!
\file passive_scatra_algorithm.cpp

\brief Transport of passive scalars in Navier-Stokes velocity field

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "passive_scatra_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::PassiveScaTraAlgorithm::PassiveScaTraAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn, const int disnum)
  :  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,disnum),
   samstart_(prbdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_START")),
   samstop_(prbdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_STOP"))
{
 // no stuff to add here at the moment
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::PassiveScaTraAlgorithm::~PassiveScaTraAlgorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::TimeLoop()
{
  // write out inital state
  // Output();

  // time loop (no-subcycling at the moment)
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // solve nonlinear Navier-Stokes system
    DoFluidStep();

    // solve transport (convection-diffusion) equations for passive scalar
    DoTransportStep();

    // update all field solvers
    Update();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n   TIME STEP     \n******************\n";
    cout<<"\nStep:   " << Step() << " / " << NStep() << "\n";
    cout<<"\n******************\n   FLUID SOLVER  \n******************\n";
  }

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
void SCATRA::PassiveScaTraAlgorithm::DoFluidStep()
{
  // solve nonlinear Navier-Stokes system
  if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    FluidField().MultiCorrector();
  else
    FluidField().NonlinearSolve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::DoTransportStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n TRANSPORT SOLVER \n******************\n";
  }

  // transfer velocities to scalar transport field solver
  // NOTE: so far, the convective velocity is chosen to equal the fluid velocity
  //       since it is not yet clear how the grid velocity should be interpolated
  //       properly -> hence, PassiveScaTraAlgorithm does not support moving
  //       meshes yet

  if (FluidField().TimIntScheme()== INPAR::FLUID::timeint_gen_alpha)
  {
    ScaTraField().SetVelocityField(
        FluidField().Velaf(),
        FluidField().Accam(),
        FluidField().Velaf(),
        Teuchos::null, // no fsvel in gammis gen alpha fluid code
        Teuchos::null,
        FluidField().Discretization());
  }
  else if (FluidField().TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
  {
    ScaTraField().SetVelocityField(
      FluidField().Velaf(),
        FluidField().Accam(),
        FluidField().Velaf(),
        FluidField().FsVel(),
        Teuchos::null,
        FluidField().Discretization());;
  }
  else
  {
    ScaTraField().SetVelocityField(
      FluidField().Velnp(),
        FluidField().Hist(),
        FluidField().Velnp(),
        FluidField().FsVel(),
        Teuchos::null,
        FluidField().Discretization()
    );
  }

  // solve the linear convection-diffusion equation(s)
  ScaTraField().Solve();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.

  if ((Step()>=samstart_) and (Step()<=samstop_))
  {
    // if statistics for one-way coupled problems is performed, provide
    // the field for the first scalar!
    FluidField().SetTimeLomaFields(
        ScaTraField().Phinp(),
        0.0,
        ScaTraField().TrueResidual(),
        ScaTraField().Discretization(),
        0 // do statistics for FIRST dof at every node!!
    );
  }

  FluidField().StatisticsAndOutput();
  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
