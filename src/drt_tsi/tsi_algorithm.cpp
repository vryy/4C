/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.cpp

\brief  Basis of all TSI algorithms that perform a coupling between the linear
        momentum equation and the heat conduction equation

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_algorithm.H"
#include "../drt_inpar/inpar_tsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(
  Epetra_Comm& comm
  )
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    tempincnp_(rcp(new Epetra_Vector(*(ThermoField().Tempnp()))))
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();

  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // decide if one-way coupling or full coupling
  INPAR::TSI::PartitionedCouplingMethod method =
    Teuchos::getIntegralValue<INPAR::TSI::PartitionedCouplingMethod>(tsidyn,"PARTITIONED");

  if (method == INPAR::TSI::OneWay)
  {
    // what kind of one-way coupling should be applied
    displacementcoupling_
      = tsidyn.get<std::string>("COUPVARIABLE") == "Displacement";

    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // build a proxy of the structure discretization for the temperature field
      Teuchos::RCP<DRT::DofSet> structdofset
        = StructureField().Discretization()->GetDofSetProxy();
      // check if ThermoField has 2 discretizations, so that coupling is possible
      if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
        dserror("unexpected dof sets in thermo field");

      // build matrices again and now consider displacements in the thermo field
      ThermoField().TSIMatrix();
    }

    else // temperature as coupling variable (deformation due to temperature change)
    {
      // build a proxy of the temperature discretization for the structure field
      Teuchos::RCP<DRT::DofSet> thermodofset
        = ThermoField().Discretization()->GetDofSetProxy();
      // check if StructField has 2 discretizations, so that coupling is possible
      if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
        dserror("unexpected dof sets in structure field");

      // now build the matrices again and consider the temperature dependence
      StructureField().TSIMatrix();
    }
  }
  // 30.07.10
  else if (method == INPAR::TSI::SequStagg || method == INPAR::TSI::IterStagg)
  {
    // build a proxy of the structure discretization for the temperature field
    Teuchos::RCP<DRT::DofSet> structdofset
      = StructureField().Discretization()->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSet> thermodofset
      = ThermoField().Discretization()->GetDofSetProxy();

    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in thermo field");
    if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
      dserror("unexpected dof sets in structure field");

    // now build the matrices again and consider dependencies to 2nd field
    ThermoField().TSIMatrix();
    StructureField().TSIMatrix();
  }

//    // now check if the two dofmaps are available and then bye bye
//    cout << "structure dofmap" << endl;
//    cout << *StructureField().DofRowMap(0) << endl;
//    cout << "thermo dofmap" << endl;
//    cout << *StructureField().DofRowMap(1) << endl;
//    cout << "thermo dofmap" << endl;
//    cout << *ThermoField().DofRowMap(0) << endl;
//    cout << "structure dofmap" << endl;
//    cout << *ThermoField().DofRowMap(1) << endl;
//    exit(0);

}


/*----------------------------------------------------------------------*
 | destructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ReadRestart(int step)
{
//  ThermoField().ReadRestart(step);
//  StructureField().ReadRestart(step);
//  SetTimeStep(ThermoField().GetTime(),step);

  return;
}


/*----------------------------------------------------------------------*
 | Timeloop                                                  dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoop()
{
  // get an idea of the temperatures (like in partitioned FSI)
  tempn_ = ThermoField().ExtractTempn();
  dispn_ = StructureField().ExtractDispn();

  // ==================================================================

  // time loop
  while (NotFinished())
  {
    // get active nodes from structural contact simulation
    RCP<MORTAR::ManagerBase> cmtman = StructureField().ContactManager();

    // tsi with or without contact
    // only tsi
    if (cmtman == Teuchos::null)
    {
      // call the TSI parameter list
      const Teuchos::ParameterList& tsidyn
        = DRT::Problem::Instance()->TSIDynamicParams();

      // decide if apply one-way coupling or full coupling
      INPAR::TSI::PartitionedCouplingMethod method =
        Teuchos::getIntegralValue<INPAR::TSI::PartitionedCouplingMethod>(tsidyn,"PARTITIONED");

      // switch TSI algorithm (synchronous)
      // only the coupling field knows the 2nd (coupling) field
      if (method == INPAR::TSI::OneWay)
        TimeLoopOneWay();
      // complete volume coupling system
      // sequential staggered scheme
//      else if (method == INPAR::TSI::SequStagg)
//        TimeLoopSequStagg();
      // iterative staggered scheme
      else if (method == INPAR::TSI::IterStagg)
        TimeLoopFull();

    } // tsi

    // thermo-structure interaction with contact in seperate routine
    else
    {
      //****************************************************************//
      // this algorithm consists of two parts within one time step      //
      // 1) solution of nonlinear structural system without influence   //
      //    of temperature.                                             //
      // 2) solution of temperature problem                             //
      //    with heat transfer over the contacting surface (as a result //
      //    from the structural problem).                               //
      //******************************************************************

      // only for frictional contact so far
      // as information are needed from the friction node
      if((cmtman->GetStrategy().Friction())==false)
        dserror ("Thermo-Structure interaction only for frictional contact so far");

      // counter and print header
      IncrementTimeAndStep();
      PrintHeader();
      
      // evaluate reference state
      StructureField().EvaluateReferenceState();

      // predict ans solve structural system
      StructureField().PrepareTimeStep();
      StructureField().Solve();
      
      // initialize contact manager of thermo field
      ThermoField().SetStructContact(cmtman,StructureField().Discretization());

      // predict and solve thermal problem
      ThermoField().PrepareTimeStep();
      ThermoField().Solve();
      
      // update all single field solvers
      Update();

      // write output to screen and files
      Output();
    }
  } // time loop

  // ==================================================================

}


/*----------------------------------------------------------------------*
 | One-way coupling (only one field knows his counterpart)   dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopOneWay()
{
    // counter and print header
    IncrementTimeAndStep();
    PrintHeader();

    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // predict the structure field without influence of temperatures
      StructureField().PrepareTimeStep();

      // Begin Nonlinear Solver / Outer Iteration ******************************

      // iterate between the two fields
      int  itnum = 0;
      bool stopnonliniter = false;

      // Outer Iteration loop starts
      if (Comm().MyPID()==0)
      {
        cout<<"\n";
        cout<<"**************************************************************\n";
        cout<<"      OUTER ITERATION LOOP \n";
        printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
          ThermoField().GetTimeNumStep());
        cout<<"**************************************************************\n";
      }

      // Start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // store temperature from first solution for convergence check
        tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

        // get current displacement due to Solve Structure Step
        const Teuchos::RCP<Epetra_Vector> dispnp = DoStructureStep();
        // 04.08.10 comment ExtractVelnp() to exclude oszillation
        // extract the velocities of the current solution
        const Teuchos::RCP<Epetra_Vector> velnp = StructureField().ExtractVelnp();
        // the displacement -> velocity conversion
//        const Teuchos::RCP<Epetra_Vector> velnp = CalcVelocity(dispnp);

        // solve thermo system
        // and therefore use the u-solution calculated in DoStructureStep
        // including: - ApplyDisplacements()
        //            - PrepareTimeStep()
        //            - Solve()
        DoThermoStep(dispnp, velnp);

        // check convergence of temperature field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

      } // end OUTER ITERATION

      // End Nonlinear Solver **************************************

      // update all single field solvers
      Update();

      // extract final displacements,
      // since we did update, this is very easy to extract
      dispn_ = StructureField().ExtractDispn();

    } // displacement coupling

    else // temperature coupling (deformation due to heating)
    {
      // predict the thermal field without influence of structure
      ThermoField().PrepareTimeStep();

      // Begin Nonlinear Solver / Outer Iteration ******************************

      // iterate between the two fields
      int  itnum = 0;
      bool stopnonliniter = false;

      // Outer Iteration loop starts
      if (Comm().MyPID()==0)
      {
        cout<<"\n";
        cout<<"**************************************************************\n";
        cout<<"      OUTER ITERATION LOOP \n";
        printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
          ThermoField().GetTimeNumStep());
        cout<<"**************************************************************\n";
      }

      // Start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // store temperature from first solution for convergence check
        tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

        // get current temperatures due to Solve ThermoStep
        const Teuchos::RCP<Epetra_Vector> tempnp = DoThermoStep();

        // solve structure system
        // and therefore use the temperature solution calculated in DoThermoStep
        // including: - ApplyTemperatures()
        //            - PrepareTimeStep()
        //            - Solve()
        DoStructureStep(tempnp);
        
        // check convergence of temperature field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

      } // end OUTER ITERATION

      // End Nonlinear Solver **************************************

      // update all single field solvers
      Update();

      // extract final temperatures,
      // since we did update, this is very easy to extract
      tempn_ = ThermoField().ExtractTempn();

    }  // temperature coupling////

      // write output to screen and files
      Output();
}


/*----------------------------------------------------------------------*
 | Full coupling (iterative staggered scheme)                dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopFull()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  // Outer iteration loop
  OuterIterationLoop();

  // update all single field solvers
  Update();

  // extract final displacement and velocity
  // since we did update, this is very easy to extract
  dispn_ = StructureField().ExtractDispn();
  tempn_ = ThermoField().ExtractTempn();

  // write output to screen and files
  Output();

}


/*----------------------------------------------------------------------*
 | Outer Loop with convergence check                         dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::OuterIterationLoop()
{
  // Begin Nonlinear Solver / Outer Iteration ******************************

  // iterate between the two fields
  int  itnum = 0;
  bool stopnonliniter = false;

  // Outer Iteration loop starts
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"**************************************************************\n";
    cout<<"      OUTER ITERATION LOOP \n";
    printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
      ThermoField().GetTimeNumStep());
    cout<<"**************************************************************\n";
  }

  // Start OUTER ITERATION
  while (stopnonliniter == false)
  {
    itnum ++;

    // store temperature from first solution for convergence check
    tempincnp_->Update(1.0,*ThermoField().Tempn(),0.0);

    // get current temperatures due to Solve ThermoStep, like predictor in FSI
    if (itnum == 1)
    {
      tempn_ = ThermoField().ExtractTempn();
    }
    else
      tempn_ = ThermoField().ExtractTempnp();

    // solve structure system and extract current displacements
    // and therefore use the temperature solution calculated in DoThermoStep
    // including: - ApplyTemperatures()
    //            - PrepareTimeStep()
    //            - Solve()
    const Teuchos::RCP<Epetra_Vector> dispnp = DoFullStructureStep(tempn_);
    // extract the velocities of the current solution
    const Teuchos::RCP<Epetra_Vector> velnp = StructureField().ExtractVelnp();

    // solve thermo system
    // and therefore use the u-solution calculated in DoStructureStep
    // including: - ApplyDisplacements()
    //            - PrepareTimeStep()
    //            - Solve()
    const Teuchos::RCP<Epetra_Vector> tempn_ = DoFullThermoStep(dispnp, velnp);

    // check convergence of temperature field for "partitioned scheme"
    stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

  } // end OUTER ITERATION
  // End Nonlinear Solver **************************************

  return;
}  // OuterIterationLoop()


/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 03/10 |
 | apply temperatures, coupling variable in structure field             |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoStructureStep(Teuchos::RCP<Epetra_Vector> temp)
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }
  // call the current temperatures
  StructureField().ApplyTemperatures(temp);

  // call the predictor here, because the temperature is considered too
  StructureField().PrepareTimeStep();

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  return;
}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 05/10 |
 | extract end temperatures for coupling to displacement field          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::DoThermoStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }
  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();
  // now extract the current temperatures and pass it to the structure
  return ThermoField().ExtractTempnp();
}


/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 05/10 |
 | extract current displacements needed for coupling to thermo field    |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>  TSI::Algorithm::DoStructureStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  // now extract the current temperatures and pass it to the structure
  return StructureField().ExtractDispnp();
}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 05/10 |
 | coupling of displacements to thermo field                            |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoThermoStep(
  Teuchos::RCP<Epetra_Vector> disp,
  Teuchos::RCP<Epetra_Vector> vel
  )
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }

  // call the current displacements and velocities
  ThermoField().ApplyStructVariables(disp,vel);

  // call the predictor here, because displacement field is considered, too
  ThermoField().PrepareTimeStep();

  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();

  return;
}

/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 06/10 |
 | extract current displacements needed for coupling to thermo field    |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>  TSI::Algorithm::DoFullStructureStep(
  Teuchos::RCP<Epetra_Vector> temp
  )
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }
  // call the current temperatures
  StructureField().ApplyTemperatures(temp);

  // call the predictor here, because the temperature is considered too
  StructureField().PrepareTimeStep();
  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();

  // now extract the current temperatures and pass it to the structure
  return StructureField().ExtractDispnp();

}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 06/10 |
 | coupling of displacements to thermo field and extract current        |
 | displacements needed for coupling to thermo field                    |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>  TSI::Algorithm::DoFullThermoStep(
  Teuchos::RCP<Epetra_Vector> disp,
  Teuchos::RCP<Epetra_Vector> vel
  )
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }

  // call the current displacements and velocities
  ThermoField().ApplyStructVariables(disp,vel);

  // call the predictor here, because displacement field is considered, too
  ThermoField().PrepareTimeStep();

  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();

  // now extract the current temperatures and pass it to the structure
  return ThermoField().ExtractTempnp();
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField().Update();
  ThermoField().Update();
  return;
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField().Output();

  // write the thermo output (temperatures at the moment) to the the structure output
  // get disc writer from structure field
  Teuchos::RCP<IO::DiscretizationWriter> output = StructureField().DiscWriter();

  // get the temperature and the noderowmap of thermo discretization
  Epetra_Vector temperature = *(ThermoField().Tempn());
  const Epetra_Map* temprowmap = ThermoField().Discretization()->NodeRowMap();

  // replace map and write it to output
  temperature.ReplaceMap(*temprowmap);
  RCP<Epetra_Vector> temp = rcp(new Epetra_Vector(temperature));
  output->WriteVector("temperature",temp);

  ThermoField().Output();
}


/*----------------------------------------------------------------------*
 | convergence check for the temperature field               dano 02/10 |
 | originally by vg 01/09                                               |
 *----------------------------------------------------------------------*/
bool TSI::Algorithm::ConvergenceCheck(
  int itnum,
  const int itmax,
  const double ittol
  )
{
  // convergence check based on the temperature increment
  bool stopnonliniter = false;

  //    | temperature increment |_2
  //  -------------------------------- < Tolerance
  //     | temperature_n+1 |_2

  // Variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double tempincnorm_L2(0.0); // pot
  double tempnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  tempincnp_->Update(1.0,*(ThermoField().Tempnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  ThermoField().Tempnp()->Norm2(&tempnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    cout<<"\n";
    cout<<"**************************************************************\n";
    cout<<"    OUTER ITERATION STEP    \n";
    cout<<"**************************************************************\n";
    printf("+--------------+------------------------+--------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]     -|--  temp-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         |",
         itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+\n");
  }

  // Converged
  if (tempincnorm_L2/tempnorm_L2 <= ittol)
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !  |\n", itnum,itmax);
      printf("+--------------+------------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax) and (tempincnorm_L2/tempnorm_L2 > ittol))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf("|     >>>>>> not converged in itemax steps!                  |\n");
      printf("+--------------+------------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Algorithm::ConvergenceCheck


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 05/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann               |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> dispnp
  ) const
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = rcp(new Epetra_Vector(*dispn_));

  // calculate velocity with timestep Dt()
  //  V_n+1 = (D_n+1 - D_n) / Dt
  vel->Update(1./Dt(), *dispnp, -1./Dt());

  return vel;
}

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
