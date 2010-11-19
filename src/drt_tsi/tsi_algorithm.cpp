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
TSI::Algorithm::Algorithm(Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    tempincnp_(rcp(new Epetra_Vector(*(ThermoField().Tempnp())))),
    dispincnp_(rcp(new Epetra_Vector(*(StructureField().Dispnp()))))
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();

  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // decide if one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields method =
    Teuchos::getIntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");
  // coupling variable
  displacementcoupling_
    = tsidyn.get<std::string>("COUPVARIABLE") == "Displacement";
  if (displacementcoupling_)
    cout << "Coupling variable: displacement" << endl;
  else
    cout << "Coupling variable: temperature" << endl;

  // if structure field is quasi-static --> CalcVelocity
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // major switch to different time integrators
  quasistatic_
    = Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP")==INPAR::STR::dyna_statics;

  if (method == INPAR::TSI::OneWay)
  {
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
 | destructor (public)                                       dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ReadRestart(int step)
{
  ThermoField().ReadRestart(step);
  StructureField().ReadRestart(step);
  SetTimeStep(ThermoField().GetTime(),step);

  return;
}


/*----------------------------------------------------------------------*
 | time loop                                                 dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoop()
{
  // get an idea of the temperatures (like in partitioned FSI)
  tempn_ = ThermoField().ExtractTempn();
  // get an idea of the displacements and velocities
  dispn_ = StructureField().ExtractDispn();
  veln_ = StructureField().ExtractVeln();

  // ==================================================================
  
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();

  // decide if apply one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields method =
    Teuchos::getIntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");
  
  // get active nodes from structural contact simulation
  RCP<MORTAR::ManagerBase> cmtman = StructureField().ContactManager();

  // tsi with or without contact
  // only tsi
  if (cmtman == Teuchos::null)
  {
    // switch TSI algorithm (synchronous)
    // only the coupling field knows the 2nd (coupling) field
    if (method == INPAR::TSI::OneWay)
      TimeLoopOneWay();
    // complete volume coupling system
    // sequential staggered scheme
    else if (method == INPAR::TSI::SequStagg)
      TimeLoopSequStagg();
    // iterative staggered scheme
    else if (method == INPAR::TSI::IterStagg)
      TimeLoopFull();
  } // tsi

  // thermo-structure interaction with contact in seperate routine
  else
  {
    // time loop
    while (NotFinished())
    {
      // only for frictional contact so far
      // as information are needed from the friction node
      if((cmtman->GetStrategy().Friction())==false)
        dserror ("Thermo-Structure interaction only for frictional contact so far");

      // counter and print header
      IncrementTimeAndStep();
      PrintHeader();

      //****************************************************************//
      // this algorithm consists of two parts within one time step      //
      // 1) solution of nonlinear structural system without influence   //
      //    of temperature.                                             //
      // 2) solution of temperature problem                             //
      //    with heat transfer over the contacting surface (as a result //
      //    from the structural problem).                               //
      //******************************************************************
      if(method == INPAR::TSI::OneWay)
      {
        // predict and solve structural system
        StructureField().PrepareTimeStep();
        StructureField().Solve();
  
        // initialize contact manager of thermo field
        ThermoField().SetStructContact(cmtman,StructureField().Discretization());
  
        ThermoField().PrepareTimeStep();
        ThermoField().Solve();
      }
      //****************************************************************//
      // this algorithm consists of the iteration between the           //   
      // two single fields until convergence is achieved. The two fields// 
      // influence each other:                                          //   
      // 1) The thermal field is influenced by the structural field with// 
      //    the contact surface (mortar matrices, active nodes....) and //
      //    frictional heating.                                         // 
      // 2) The structural field is influenced by the thermal field     //   
      //    with a temperature dependent material law.                  //
      //******************************************************************
      else if (method == INPAR::TSI::IterStagg)
      {
        // prepare time step
        ThermoField().PrepareTimeStep();
        StructureField().PrepareTimeStep();
        
        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;
        
        while (stopnonliniter == false)
        {
          itnum ++;
              
          // store temperature from first solution for convergence check (like in
          // elch_algorithm: use current values)
          tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
          dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);
        
          // initialize contact manager of thermo field
          ThermoField().SetStructContact(cmtman,StructureField().Discretization());
        
          ThermoField().PrepareTimeStep();
          ThermoField().Solve();

          // extract current temperature field
          const Teuchos::RCP<Epetra_Vector> tempnp = ThermoField().ExtractTempnp();
          
          // apply current temperatur field
          StructureField().ApplyTemperatures(tempnp);
          
          // solve structure system
          StructureField().Solve();

          // check convergence of temperature field for "partitioned scheme"
          stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);
        }
      }
      else
        dserror("No sequential staggered coupling algorithm with contact");
 
      // update all single field solvers
      Update();

      // write output to screen and files
      Output();
    } // time loop
  }

  // ==================================================================

}


/*----------------------------------------------------------------------*
 | One-way coupling (only one field knows his counterpart)   dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopOneWay()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    IncrementTimeAndStep();
    PrintHeader();

    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // predict the structure field without influence of temperatures
      StructureField().PrepareTimeStep();

      // get current displacement due to Solve Structure Step
      const Teuchos::RCP<Epetra_Vector> dispnp = DoStructureStep();

      // extract the velocities of the current solution
      Teuchos::RCP<Epetra_Vector> velnp;
      // major switch to different time integrators
      if (quasistatic_)
      {
        // the displacement -> velocity conversion
        // Use this kind of calculation for quasi-static structure calculation
        velnp = CalcVelocity(dispnp);
      }
      else
      {
        velnp = StructureField().ExtractVelnp();
      }

      // solve thermo system
      // and therefore use the u-solution calculated in DoStructureStep
      // including: - ApplyDisplacements()
      //            - PrepareTimeStep()
      //            - Solve()
      DoThermoStep(dispnp, velnp);

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

      // get current temperatures due to Solve ThermoStep
      const Teuchos::RCP<Epetra_Vector> tempnp = DoThermoStep();

      // solve structure system
      // and therefore use the temperature solution calculated in DoThermoStep
      // including: - ApplyTemperatures()
      //            - PrepareTimeStep()
      //            - Solve()
      DoStructureStep(tempnp);

      // update all single field solvers
      Update();

      // extract final temperatures,
      // since we did update, this is very easy to extract
      tempn_ = ThermoField().ExtractTempn();

    }  // temperature coupling

    // write output to screen and files
    Output();

  } // time loop
}

/*----------------------------------------------------------------------*
 | One-way coupling between the fields                       dano 09/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopSequStagg()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    IncrementTimeAndStep();
    PrintHeader();

    // like implicit-implicit staggered scheme, compare Farhat & Park, 1991
    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // Begin Nonlinear Solver / Outer Iteration ******************************

      // Predict the displacement field
      // get structure variables of old time step (d_n, v_n)
      dispn_ = StructureField().ExtractDispn();
      if(Step()==0) veln_ = StructureField().ExtractVeln();

      // Solve thermo field with predicted structural values (Apply d_n,v_n to
      // THR, Solve coupled equation, extract new temperatures T_n+1)
      tempn_ = DoFullThermoStep(dispn_, veln_);

      // Correct structure field with just calculates temperatures T_n+1 (Apply
      // T_n+1 to STR, Solve coupled equation
      const Teuchos::RCP<Epetra_Vector> dispnp = DoFullStructureStep(tempn_);

      // End Nonlinear Solver **************************************

      if(quasistatic_)
        veln_ = CalcVelocity(StructureField().ExtractDispnp());
      else
        veln_ = StructureField().ExtractVelnp();

      // update all single field solvers
      Update();

      // extract final displacements,
      // since we did update, this is very easy to extract
      dispn_ = StructureField().ExtractDispn();

    } // end displacement coupling

    // temperature is coupling variable
    else  // (deformation due to heating)
    {
      // Begin Nonlinear Solver / Outer Iteration ******************************

      // get temperature solution of old time step
      tempn_ = ThermoField().ExtractTempn();

      // get current displacement due to Solve Structure Step
      const Teuchos::RCP<Epetra_Vector> dispnp = DoFullStructureStep(tempn_);
      // extract the velocities of the current solution
      // the displacement -> velocity conversion
      if(quasistatic_)
        velnp_ = CalcVelocity(dispnp);
      // quasistatic to exlude oszillation
      else
        velnp_ = StructureField().ExtractVelnp();

      const Teuchos::RCP<Epetra_Vector> tempnp = DoFullThermoStep(dispnp, velnp_);

      // End Nonlinear Solver **************************************

      // update all single field solvers
      Update();

      // extract final displacements,
      // since we did update, this is very easy to extract
      tempn_ = ThermoField().ExtractTempn();
    } // temperature coupling

    // write output to screen and files
    Output();

  } // time loop

} // End TimeLoopSequStagg


/*----------------------------------------------------------------------*
 | Full coupling (iterative staggered scheme)                dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopFull()
{
  // time loop
  while (NotFinished())
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

  } // time loop

} // TimeLoopFull


/*----------------------------------------------------------------------*
 | Outer Loop with convergence check                         dano 10/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::OuterIterationLoop()
{
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

  // displacement predictor
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    // mechanical predictor for the coupling iteration outside the loop
    // get structure variables of old time step (d_n, v_n)
    // d^p_n+1 = d_n, v^p_n+1 = v_n
    dispn_ = StructureField().ExtractDispn();
    Teuchos::RCP<Epetra_Vector> dispnp;
    if (Step()==1)
    {
      dispnp = StructureField().ExtractDispn();
      veln_ = StructureField().ExtractVeln();
    }
    // else the velocity of the last converged step

    // Start OUTER ITERATION
    while (stopnonliniter == false)
    {
      itnum ++;

      // store temperature from first solution for convergence check (like in
      // elch_algorithm: use current values)
      tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
      dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);

      // kind of mechanical predictor
      // 1st iteration: get structure variables of old time step (d_n, v_n)
      if (itnum == 1) dispnp = StructureField().ExtractDispn();
      // set dispnp for the next time step
      else if (itnum != 1) dispnp = StructureField().ExtractDispnp();

      // Begin Nonlinear Solver / Outer Iteration ******************************

      // Solve thermo field with predicted structural values (Apply d_n,v_n to
      // THR, Solve coupled equation, extract new temperatures T_n+1)
      tempn_ = DoFullThermoStep(dispnp, veln_);

      // Correct structure field with just calculates temperatures T_n+1 (Apply
      // T_n+1 to STR, Solve coupled equation
      dispnp = DoFullStructureStep(tempn_);

      // End Nonlinear Solver **************************************

      if(quasistatic_)
        veln_ = CalcVelocity(dispnp);
      else
        veln_ = StructureField().ExtractVelnp();

      // check convergence of both field for "partitioned scheme"
      stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

    } // end OUTER ITERATION
  }  // end displacement predictor

  // temperature is predictor
  else  // (deformation due to heating)
  {
    // thermal predictor for the coupling iteration outside the loop
    // get temperature of old time step (T_n)
    // T^p_n+1 = T_n
    tempn_ = ThermoField().ExtractTempn();

    // Start OUTER ITERATION
    while (stopnonliniter == false)
    {
      itnum ++;
      // get current temperatures due to Solve ThermoStep, like predictor in FSI
      if (itnum == 1){ tempn_ = ThermoField().ExtractTempn(); }
      else { tempn_ = ThermoField().ExtractTempnp(); }

      // store temperature from first solution for convergence check (like in
      // elch_algorithm: here the current values are needed)
      tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
      dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);

      // solve structure system and extract current displacements
      // and therefore use the temperature solution calculated in DoThermoStep
      // including: - ApplyTemperatures()
      //            - PrepareTimeStep()
      //            - Solve()
      const Teuchos::RCP<Epetra_Vector> dispnp = DoFullStructureStep(tempn_);
      // extract the velocities of the current solution
      // the displacement -> velocity conversion
      if(quasistatic_)
        velnp_ = CalcVelocity(dispnp);
      // quasistatic to exlude oszillation
      else
        velnp_ = StructureField().ExtractVelnp();

      // solve thermo system
      // and therefore use the u-solution calculated in DoStructureStep
      // including: - ApplyDisplacements()
      //            - PrepareTimeStep()
      //            - Solve()
      tempn_ = DoFullThermoStep(dispnp, velnp_);

      // End Nonlinear Solver **************************************

      // check convergence of both field for "partitioned scheme"
      stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

    } // end OUTER ITERATION

  } // temperature predictor

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
 | convergence check for both fields (thermo & structure)    dano 02/10 |
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
  double tempincnorm_L2(0.0);
  double tempnorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  tempincnp_->Update(1.0,*(ThermoField().Tempnp()),-1.0);
  dispincnp_->Update(1.0,*(StructureField().Dispnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  ThermoField().Tempnp()->Norm2(&tempnorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  StructureField().Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    cout<<"\n";
    cout<<"***********************************************************************************\n";
    cout<<"    OUTER ITERATION STEP    \n";
    cout<<"***********************************************************************************\n";
    printf("+--------------+------------------------+--------------------+--------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]     -|--  temp-inc      --|--  disp-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |",
         itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+--------------------+\n");
  }

  // Converged
  if ((tempincnorm_L2/tempnorm_L2 <= ittol) &&
      (dispincnorm_L2/dispnorm_L2 <= ittol))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax);
      printf("+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax) and
       ((tempincnorm_L2/tempnorm_L2 > ittol) || (dispincnorm_L2/dispnorm_L2 > ittol))
     )
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf("|     >>>>>> not converged in itemax steps!                                       |\n");
      printf("+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Algorithm::ConvergenceCheck


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 05/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
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
