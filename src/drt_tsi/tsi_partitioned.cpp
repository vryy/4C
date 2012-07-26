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
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_partitioned.H"
#include "tsi_defines.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_adapter/adapter_thermo.H"

// contact
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_mortar/mortar_manager_base.H"

//! Note: The order of calling the two BasePartitioned-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Partitioned::Partitioned(const Epetra_Comm& comm)
: Algorithm(comm),
  tempincnp_(rcp(new Epetra_Vector(*(ThermoField()->Tempnp())))),
  dispincnp_(rcp(new Epetra_Vector(*(StructureField()->Dispnp())))),
  disp_(Teuchos::null),
  veln_(Teuchos::null),
  velnp_(Teuchos::null),
  temp_(rcp(new Epetra_Vector(*(ThermoField()->Tempn())))),
  del_(Teuchos::null),
  delhist_(Teuchos::null),
  omegan_(0.0),
  omeganp_(0.0)
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();

  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // decide if one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields coupling =
    DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");
  // coupling variable
  displacementcoupling_
    = tsidyn.get<std::string>("COUPVARIABLE")=="Displacement";
  if (displacementcoupling_)
    cout << "Coupling variable: displacement" << endl;
  else
    cout << "Coupling variable: temperature" << endl;

  // if structure field is quasi-static --> CalcVelocity
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // major switch to different time integrators
  quasistatic_
    = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP")==INPAR::STR::dyna_statics;

  // choose algorithm depending on solution type
  switch (coupling)
  {
  case INPAR::TSI::OneWay :
  {
    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // build a proxy of the structure discretization for the temperature field
      Teuchos::RCP<DRT::DofSet> structdofset
        = StructureField()->Discretization()->GetDofSetProxy();
      // check if ThermoField has 2 discretizations, so that coupling is possible
      if (ThermoField()->Discretization()->AddDofSet(structdofset)!=1)
        dserror("unexpected dof sets in thermo field");
    }

    else // temperature as coupling variable (deformation due to temperature change)
    {
      // build a proxy of the temperature discretization for the structure field
      Teuchos::RCP<DRT::DofSet> thermodofset
        = ThermoField()->Discretization()->GetDofSetProxy();
      // check if StructField has 2 discretizations, so that coupling is possible
      if (StructureField()->Discretization()->AddDofSet(thermodofset)!=1)
        dserror("unexpected dof sets in structure field");
    }
    break;
  }
  case INPAR::TSI::SequStagg :
  case INPAR::TSI::IterStagg :
  case INPAR::TSI::IterStaggAitken :
  case INPAR::TSI::IterStaggAitkenIrons :
  {
    // build a proxy of the structure discretization for the temperature field
    Teuchos::RCP<DRT::DofSet> structdofset
      = StructureField()->Discretization()->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSet> thermodofset
      = ThermoField()->Discretization()->GetDofSetProxy();

    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (ThermoField()->Discretization()->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in thermo field");
    if (StructureField()->Discretization()->AddDofSet(thermodofset)!=1)
      dserror("unexpected dof sets in structure field");
    break;
  }
  default:
     dserror("Unknown solutiontype for thermo-structure interaction: %d",coupling);
  }  // end switch

#ifdef TSIPARTITIONEDASOUTPUT
    // now check if the two dofmaps are available and then bye bye
    cout << "structure dofmap" << endl;
    cout << *StructureField()->DofRowMap(0) << endl;
    cout << "thermo dofmap" << endl;
    cout << *StructureField()->DofRowMap(1) << endl;
    cout << "thermo dofmap" << endl;
    cout << *ThermoField()->DofRowMap(0) << endl;
    cout << "structure dofmap" << endl;
    cout << *ThermoField()->DofRowMap(1) << endl;
//    exit(0);
#endif // TSIPARTITIONEDASOUTPUT

    // contact
    if(StructureField()->ContactManager() != null)
      ThermoField()->PrepareThermoContact(StructureField()->ContactManager(),StructureField()->Discretization());

}  // Constructor


/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Partitioned::~Partitioned()
{
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::ReadRestart(int step)
{
  ThermoField()->ReadRestart(step);
  StructureField()->ReadRestart(step);
  SetTimeStep(ThermoField()->GetTime(),step);

  return;
}  // ReadRestart


/*----------------------------------------------------------------------*
 | prepare time step (public)                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::PrepareTimeStep()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

}  // PrepareTimeStep()


/*----------------------------------------------------------------------*
 | initialise internal variables needed as guess             dano 02/12 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::SetupSystem()
{
  // get an idea of the temperatures (like in partitioned FSI)
  temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
  // get an idea of the displacements and velocities
  disp_ = StructureField()->ExtractDispn();
  veln_ = StructureField()->ExtractVeln();

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | time loop                                                 dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::TimeLoop()
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();

  // decide if apply one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields coupling =
    DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // get active nodes from structural contact simulation
  Teuchos::RCP<MORTAR::ManagerBase> cmtman = StructureField()->ContactManager();

  // tsi with or without contact
  // only tsi
  if (cmtman == Teuchos::null)
  {
    // time loop
    while (NotFinished())
    {
      // counter and print header
      PrepareTimeStep();

      // normally PrepareTimeStep() should be called outside the nonlinear loop,
      // but at this point the coupling variables are not known, yet!
      // call PrepareTimeStep() after ApplyTemperatures()/ApplyStructVariables()

      // switch TSI algorithm (synchronous)
      // only the one field knows the 2nd field, that contains the coupling
      // variable
      // choose algorithm depending on solution type
      switch (coupling)
      {
      case INPAR::TSI::OneWay :
      {
        TimeLoopOneWay();
        break;
      }
      // complete volume coupling system
      // sequential staggered scheme
      case INPAR::TSI::SequStagg :
      {
        TimeLoopSequStagg();
        break;
      }
      // iterative staggered scheme
      case INPAR::TSI::IterStagg :
      case INPAR::TSI::IterStaggAitken :
      case INPAR::TSI::IterStaggAitkenIrons :
      {
        TimeLoopFull();
        break;
      }
      default :
        dserror("desired type of thermo-structure interaction algorithm not supported");
      }  // end switch

      // calculate stresses, strains, energies
      PrepareOutput();

      // update all single field solvers
      Update();

      // write output to screen and files
      Output();

    }  // time loop

  }  // tsi

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
      PrepareTimeStep();

      //****************************************************************//
      // this algorithm consists of two parts within one time step      //
      // 1) solution of nonlinear structural system without influence   //
      //    of temperature.                                             //
      // 2) solution of temperature problem                             //
      //    with heat transfer over the contacting surface (as a result //
      //    from the structural problem).                               //
      //******************************************************************
      if(coupling==INPAR::TSI::OneWay)
      {
        // predict and solve structural system
        StructureField()->PrepareTimeStep();
        StructureField()->Solve();

        ThermoField()->PrepareTimeStep();
        ThermoField()->Solve();
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
      else if (coupling==INPAR::TSI::IterStagg)
      {
        // prepare time step
        ThermoField()->PrepareTimeStep();
        StructureField()->PrepareTimeStep();

        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;

        while (stopnonliniter==false)
        {
          itnum ++;

          // store temperature from first solution for convergence check (like in
          // elch_algorithm: use current values)
          tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
          dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

          ThermoField()->PrepareTimeStep();
          ThermoField()->Solve();

          // extract current temperature field
          const Teuchos::RCP<Epetra_Vector> tempnp = ThermoField()->ExtractTempnp();

          // apply current temperatur field
          StructureField()->ApplyTemperatures(tempnp);

          // solve structure system
          StructureField()->Solve();

          // check convergence of temperature field for "partitioned scheme"
          stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);
        }
      }
      else
        dserror("No sequential staggered coupling algorithm with contact");

      // calculate stresses, strains, energies
      PrepareOutput();

      // update all single field solvers
      Update();

      // write output to screen and files
      Output();
    } // time loop
  }

  // ==================================================================

}  // TSI::Partitioned::TimeLoop()


/*----------------------------------------------------------------------*
 | One-way coupling (only one field knows his counterpart)   dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::TimeLoopOneWay()
{
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    /// structure field

    // predict the structure field without influence of temperatures
    StructureField()->PrepareTimeStep();

    /// solve structural system
    // get current displacement due to solve structure step
    DoStructureStep();

    // now extract the current displacements/velocities
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->ExtractDispnp();

    // extract the velocities of the current solution
    Teuchos::RCP<Epetra_Vector> velnp;
    // major switch to different time integrators
    if (quasistatic_)
    {
      // the displacement -> velocity conversion
      // use this kind of calculation for quasi-static structure calculation
      velnp = CalcVelocity(dispnp);
    }
    else
    {
      velnp = StructureField()->ExtractVelnp();
    }

    /// thermo field

    // use the structural solution calculated in DoStructureStep

    // call the current displacements and velocities and pass it to THR
    ThermoField()->ApplyStructVariables(dispnp, velnp);

    // call the predictor here, because the displacement field is now also
    // available
    ThermoField()->PrepareTimeStep();

    /// solve temperature system
    // do the solve for the time step. All boundary conditions have been set
    DoThermoStep();

    // extract final displacements,
    // since we did update, this is very easy to extract
    disp_ = StructureField()->ExtractDispnp();

  } // displacement coupling

  else // temperature coupling (deformation due to heating)
  {
    /// thermo field

    // predict the thermal field without influence of structure
    ThermoField()->PrepareTimeStep();

    /// solve temperature system
    // get current temperatures due to Solve Thermo Step
    DoThermoStep();

    // now extract the current temperatures and pass it to the structure
    const Teuchos::RCP<Epetra_Vector> tempnp = ThermoField()->ExtractTempnp();

    /// structure field

    // pass the current temperatures to the structural field
    StructureField()->ApplyTemperatures(tempnp);

    // call the predictor here, because the temperature is now also available
    StructureField()->PrepareTimeStep();

    /// solve structure system
    // do the nonlinear solve for the time step. All boundary conditions have
    // been set.
    DoStructureStep();

    // extract final temperatures,
    // since we did update, this is very easy to extract
    temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

  }  // temperature coupling

}  // TSI::Partitioned::TimeLoopOneWay()


/*----------------------------------------------------------------------*
 | One-way coupling between the fields                       dano 09/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::TimeLoopSequStagg()
{
  // like implicit-implicit staggered scheme, compare Farhat & Park, 1991
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    // begin nonlinear solver ************************************************

    // predict the displacement field (structural predictor for TSI)

    // get structure variables of old time step (d_n, v_n)
    disp_ = StructureField()->ExtractDispn();
    if( Step()==0 ) veln_ = StructureField()->ExtractVeln();

    /// thermo field

    // solve thermo field with predicted structural values (Apply d_n,v_n to
    // THR, solve coupled equation, extract new temperatures T_n+1)

    // pass the current displacements and velocities to the thermo field
    ThermoField()->ApplyStructVariables(disp_,veln_);

    // prepare time step with coupled variables
    ThermoField()->PrepareTimeStep();

    /// solve temperature system

    // do the solve for the time step. All boundary conditions have been set.
    DoThermoStep();

    // now extract the current temperatures and pass it to the structure
    temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

    /// structure field

    // pass the current temperatures to the structure field
    StructureField()->ApplyTemperatures(temp_);

    // prepare time step with coupled variables
    StructureField()->PrepareTimeStep();

    // solve coupled equation
    DoStructureStep();

    // end nonlinear solver **************************************************

    // get the velocities needed as predictor for the thermo field for the next
    // time step
    if(quasistatic_)
      veln_ = CalcVelocity(StructureField()->ExtractDispnp());
    else
      veln_ = StructureField()->ExtractVelnp();

    // extract final displacements,
    // update is called afterwards, so we extract the newest solution
    disp_ = StructureField()->ExtractDispnp();

  } // end displacement coupling

  // temperature is coupling variable
  else  // (deformation due to heating)
  {
    // begin nonlinear solver ************************************************

    // predict the temperature field (thermal predictor for TSI)

    // get temperature solution of old time step
    temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);

    /// structure field

    // get current displacement due to solve structure step

    // pass the current temperatures to the structure field
    StructureField()->ApplyTemperatures(temp_);

    // prepare time step with coupled variables
    StructureField()->PrepareTimeStep();

    // solve structural coupled equation
    DoStructureStep();

    // extract current displacements
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->ExtractDispnp();

    // extract the velocities of the current solution
    // the displacement -> velocity conversion
    if(quasistatic_)
      velnp_ = CalcVelocity(dispnp);
    // quasistatic to exlude oszillation
    else
      velnp_ = StructureField()->ExtractVelnp();

    /// thermo field

    // pass the current displacements and velocities to the thermo field
    ThermoField()->ApplyStructVariables(dispnp,velnp_);

    // prepare time step with coupled variables
    ThermoField()->PrepareTimeStep();

    /// solve temperature system
    /// do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    DoThermoStep();

    // end nonlinear solver **************************************************

    // extract final displacements,
    // update is called afterwards, so we extract the newest solution
    temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

  } // temperature coupling

} // TSI::Partitioned::TimeLoopSequStagg()


/*----------------------------------------------------------------------*
 | Full coupling (iterative staggered scheme)                dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::TimeLoopFull()
{
  // outer iteration loop
  OuterIterationLoop();

  // extract final displacement and velocity
  // update is called afterwards, so we extract the newest solution
  disp_ = StructureField()->ExtractDispnp();
  temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

} // TSI::Partitioned::TimeLoopFull()


/*----------------------------------------------------------------------*
 | Outer Loop with convergence check                         dano 10/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::OuterIterationLoop()
{
  // iterate between the two fields
  int  itnum = 0;
  bool stopnonliniter = false;

  // outer iteration loop starts
  if (Comm().MyPID()==0 and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
  {
    cout<<"\n";
    cout<<"**************************************************************\n";
    cout<<"      OUTER ITERATION LOOP \n";
    printf("      Time Step %3d/%3d \n",ThermoField()->GetTimeStep(),
      ThermoField()->GetTimeNumStep());
    cout<<"**************************************************************\n";
  }

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();
  // decide if one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields coupling =
    DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // Pure iterative staggered algorithms
  // iterative staggered TSI withOUT Aitken relaxation
  if ( coupling==INPAR::TSI::IterStagg)
  {
    // structural predictor
    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // mechanical predictor for the coupling iteration outside the loop
      // get structure variables of old time step (d_n, v_n)
      // d^p_n+1 = d_n, v^p_n+1 = v_n
      // initialise new time step n+1 with values of old time step n
      Teuchos::RCP<Epetra_Vector> dispnp
        = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
      if ( Step()==1 )
      {
        dispnp->Update(1.0, *(StructureField()->ExtractDispn()),0.0);
        veln_ = StructureField()->ExtractVeln();
      }
      // else: use the velocity of the last converged step

      // start OUTER ITERATION
      while (stopnonliniter==false)
      {
        itnum ++;

        // store temperature from first solution for convergence check (like in
        // elch_algorithm: use current values)
        tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
        dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

        // kind of mechanical predictor
        // 1st iteration: get structure variables of old time step (d_n, v_n)
        if (itnum==1)
          dispnp->Update(1.0,*(StructureField()->ExtractDispn()),0.0);
        // for itnum>1 use the current solution dispnp of old iteration step

        // Begin Nonlinear Solver / Outer Iteration ******************************

        // thermo field

        // pass the current displacements and velocities to the thermo field
        ThermoField()->ApplyStructVariables(dispnp,veln_);

        // prepare time step with coupled variables
        if (itnum==1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        /// solve temperature system

        /// do the nonlinear solve for the time step. All boundary conditions have
        /// been set.
        DoThermoStep();

        // now extract the current temperatures and pass it to the structure
        temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

        /// structure field

        // pass the current temperatures to the structure field
        // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
        // just calculated incremental solutions
        StructureField()->ApplyTemperatures(temp_);

        // prepare time step with coupled variables
        if (itnum==1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // extract current displacements
        dispnp->Update(1.0,*(StructureField()->ExtractDispnp()),0.0);

        // end nonlinear solver / outer iteration ********************************

        if(quasistatic_)
          veln_ = CalcVelocity(dispnp);
        else
          veln_ = StructureField()->ExtractVelnp();

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
      temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);

      // start OUTER ITERATION
      while (stopnonliniter==false)
      {
        itnum ++;
        // get current temperatures due to solve thermo step, like predictor in FSI
        if (itnum==1)
        {
          temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
        }
        else
        {
          temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);
        }

        // store temperature from first solution for convergence check (like in
        // elch_algorithm: here the current values are needed)
        tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
        dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

        // begin nonlinear solver / outer iteration ******************************

        /// structure field

        // pass the current temperatures to the structure field
        // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
        // just calculated incremental solutions
        StructureField()->ApplyTemperatures(temp_);

        // prepare time step with coupled variables
        if (itnum==1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // Extract current displacements
        const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->ExtractDispnp();

        // extract the velocities of the current solution
        // the displacement -> velocity conversion
        if(quasistatic_)
          velnp_ = CalcVelocity(dispnp);
        // quasistatic to exlude oszillation
        else
          velnp_ = StructureField()->ExtractVelnp();

        /// thermo field

        // pass the current displacements and velocities to the thermo field
        ThermoField()->ApplyStructVariables(dispnp,velnp_);

        // prepare time step with coupled variables
        if (itnum==1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        /// solve coupled thermal system
        /// do the solve for the time step. All boundary conditions have been set.
        DoThermoStep();

        // end nonlinear solver / outer iteration ********************************

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

      } // end OUTER ITERATION

    } // temperature predictor

  } // iterstagg WITHOUT relaxation

  // notation according to Aitken relaxation of FSI and Mok
  // 1. calculate an Aitken factor omega
  // 2. relaxation factor relax = 1-omega
  // 3. T^{i+1} = T^i + relax^{i+1} * ( T^{i+1} - T^i )
  // 4. limit Aitken factor omega for next time step with 1.0
  else if ( coupling==INPAR::TSI::IterStaggAitken)
  {
    cout << "Iterative staggered TSI with Aitken relaxation" << endl;

    // relax the temperatures
    if (!displacementcoupling_)
    {
      if ( Step()==1 )
      {
        // thermal predictor for the coupling iteration outside the loop
        // get temperature of old time step (T_n)
        // T^p_n+1 = T_n
        temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
      }

      // initalise relaxation parameter for each iteration
      double relax = 0.0;

      // start OUTER ITERATION
      while (stopnonliniter==false)
      {
        itnum ++;

        // store temperature from first solution for convergence check (like in
        // elch_algorithm: use current values)
        // n: time indices, i: Newton iteration
        // calculate increment Delta T^{i+1}_{n+1}
        // 1. iteration step (i=0): Delta T^{n+1}_{1} = T^{n},
        //                          so far no solving has occured: T^{n+1} = T^{n}
        // i+1. iteration step:     Inc T^{i+1}_{n+1} = T^{i+1}_{n+1} - T^{i}_{n+1}
        //                     fill Inc T^{i+1}_{n+1} = T^{i}_{n+1}
        tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
        dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

        // get current temperatures due to solve thermo step, like predictor in FSI
        // 1. iteration: get temperatures of old time step (T_n)
        if (itnum==1)
        {
          temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
        }
        // else: use relaxed solution vector temp_

        // begin nonlinear solver / outer iteration ******************************

        /// structure field

        // pass the current temperatures to the structure field
        // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
        // just calculated incremental solutions
        StructureField()->ApplyTemperatures(temp_);

        // prepare time step with coupled variables
        if (itnum==1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // Extract current displacements
        const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->ExtractDispnp();

        // extract the velocities of the current solution
        // the displacement -> velocity conversion
        if(quasistatic_)
          velnp_ = CalcVelocity(dispnp);
        // quasistatic to exlude oscillation
        else
          velnp_ = StructureField()->ExtractVelnp();

        /// thermo field

        // pass the current displacements and velocities to the thermo field
        ThermoField()->ApplyStructVariables(dispnp,velnp_);

        // prepare time step with coupled variables
        if (itnum==1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        /// solve coupled thermal system
        /// do the solve for the time step. All boundary conditions have been set.
        DoThermoStep();

        // end nonlinear solver / outer iteration ********************************

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        // ------------------------------------------------------------
        // if r_{i+1} not converged in ConvergenceCheck()
        // --> apply Aitken relaxation to displacements
        // as implemented in FSI by Irons and Tuck (1969)

        // Aitken factor
        // relaxation can start when two vectors are available, starting at i=2
        // (also stated in Diss Uli, p. 82ff)
        // we need vectors from iteration step 0,1,2 to calculate relaxed step
        // with omega^i_0 = 0.0
        // Irons & Tuck:
        // increment r := old - new

        // BACI:
        // increment inc := new - old
        // omega^{i+1}_{n+1} = omega^i_{n+1} +
        //
        //                         ( r^i_{n+1} - r^{i+1}_{n+1} )^T . ( r^{i+1}_{n+1} )
        // + (omega^i_{n+1} - 1) . -------------------------------------------------
        //                                 | r^i_{n+1} - r^{i+1}_{n+1} |^2
        //
        // | r^i_{n+1} - r^{i+1}_{n+1} |^2
        // = | (-1)^2*(r^{i+1}_{n+1} - r^i_{n+1} |^2
        // = | (r^{i+1}_{n+1} - r^i_{n+1} |^2

        // initialise increment vector with solution of last iteration (i)
        // update del_ with current residual vector
        // difference of last two solutions
        if (del_ == Teuchos::null)
        {
          del_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
          delhist_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
          del_->PutScalar(1.0e20);
          delhist_->PutScalar(0.0);
        }

        // calculate difference of current (i+1) and old (i) residual vector
        // del = r^{i+1}_{n+1}
        del_->Update(1.0, *tempincnp_,0.0);
        // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
        delhist_->Update(1.0, *del_, (-1.0));
        double normdel = 0.0;
        double dot = 0.0;
        delhist_->Norm2(&normdel);
        // calculate dot product
        // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
        del_->Dot(*delhist_,&dot);

        // Aikten factor
        // omega^{i+1}_{n+1} == omeganp_
        // omega^{i}_{n+1} == omegan_
        // ome^{i+1} = ome^i + (ome^i -1) . (r^i - r^{i+1})^T . r^{i+1} / |r^{i+1} - r^{i}|^2
        omeganp_ = omegan_ + (omegan_ - 1.0) * (-dot) / (normdel*normdel);

        // relaxation parameter
        // omega_relax^{i+1} = 1- omega^{i+1}
        relax = 1.0 - omeganp_;

        // relax temperature solution for next iteration step
        // overwrite temp_ with relaxed solution vector
        // T^{i+1} = relax . T^{i+1} + (1- relax^{i+1}) T^i
        //         = T^i + relax^{i+1} * ( T^{i+1} - T^i )
        temp_->Update(relax,*del_,1.0);
        // alternatively choose a fix relaxation parameter, here e.g. 0.01
        // temp_->Update(0.01,*del_,1.0);

        // update Aitken parameter omega^{i+1}_{n+1}
        omegan_ = omeganp_;

        // update history vector with residual displacement of old iteration step
        delhist_->Update(1.0,*del_,0.0);

        // end Aitken relaxation
        // ------------------------------------------------------------

      } // end OUTER ITERATION

      // initial guess for next time step n+1
      // use maximum between omeganp_ and 1.0 as start value for omega_n+1^{i=0}
      // in case of doubt use 1.0, meaning that direction of new solution vector
      // dispnp better old one
      omegan_ = max(omeganp_, 1.0);

    } // relax temperatures
    else  // relax mechanical variables
      dserror("Relaxation of mechanical variables not yet implemented.");

  } // iterative staggered TSI with Aitken relaxation

  // another notation of relaxation according to Paper by Irony & Tuck (1969)
  // 1. calculate an Aitken factor omega == relaxation factor
  // 2. d^{i+1} = (1 - omega^{i+1}) d^{i+1} + omega^{i+1} d^i
  // 3. limit Aitken factor omega for next time step with 0.0
  else if ( coupling==INPAR::TSI::IterStaggAitkenIrons)
  {
    cout << "Iterative staggered TSI with Aitken relaxation according to Irons & Tuck" << endl;

    // relax the mechanical variables
    if (displacementcoupling_)
    {
      cout << "Displacements are relax. Be careful: not consistent, because "
        "velocities are not relaxed." << endl;

      // mechanical predictor for the coupling iteration outside the loop
      // get structure variables of old time step (d_n, v_n)
      // d^p_n+1 = d_n, v^p_n+1 = v_n
      Teuchos::RCP<Epetra_Vector> dispnp
        = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
      if ( Step()==1 )
      {
        dispnp->Update(1.0,*(StructureField()->ExtractDispn()),0.0);
        veln_ = StructureField()->ExtractVeln();
      }
      // else: use the velocity of the last converged step

      // start OUTER ITERATION
      while (stopnonliniter==false)
      {
        itnum ++;

        // store temperature from first solution for convergence check (like in
        // elch_algorithm: use current values)
        // n: time indices, i: Newton iteration
        // calculate increment Delta T^{i+1}_{n+1}
        // 1. iteration step (i=0): Delta T^{n+1}_{1} = T^{n},
        //                          so far no solving has occured: T^{n+1} = T^{n}
        // i+1. iteration step:     Inc T^{i+1}_{n+1} = T^{i+1}_{n+1} - T^{i}_{n+1}
        //                     fill Inc T^{i+1}_{n+1} = T^{i}_{n+1}
        tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
        dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

        // kind of mechanical predictor
        // 1. iteration: get structure variables of old time step (d_n, v_n)
        if (itnum==1)
          dispnp->Update(1.0,*(StructureField()->ExtractDispn()),0.0);
        // for itnum>1 use the current solution dispnp of old iteration step

        // Begin Nonlinear Solver / Outer Iteration ******************************

        // thermo field

        // pass the current displacements and velocities to the thermo field
        ThermoField()->ApplyStructVariables(dispnp,veln_);

        // prepare time step with coupled variables
        if (itnum==1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        /// solve temperature system

        /// do the nonlinear solve for the time step. All boundary conditions have
        /// been set.
        DoThermoStep();

        // now extract the current temperatures and pass it to the structure
        temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

        /// structure field

        // pass the current temperatures to the structure field
        // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
        // just calculated incremental solutions
        StructureField()->ApplyTemperatures(temp_);

        // prepare time step with coupled variables
        if (itnum==1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // check convergence of both field for "partitioned scheme"
        // calculate residual vectors
        // ~ means: not relaxed, newest solution vector after DoStructureStep()
        // r^{i+1} := Delta d^{i+1} = d^{i+1}^{~} - d^i
        // r^{i+1} == dispincnp_
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        // ------------------------------------------------------------
        // if r_{i+1} not converged in ConvergenceCheck()
        // --> apply Aitken relaxation to displacements
        // as implemented in FSI by Irons and Tuck (1969)

        // Aitken factor
        // relaxation can start when two vectors are available, starting at i=2
        // (also stated in Diss Uli, p. 82ff)
        // we need vectors from iteration step 0,1,2 to calculate relaxed step
        // with omega^i_0 = 0.0
        // Irons & Tuck:
        // increment r := old - new

        // BACI:
        // increment inc := new - old
        // omega^{i+1}_{n+1} = omega^i_{n+1} +
        //
        //                         ( r^i_{n+1} - r^{i+1}_{n+1} )^T . ( r^{i+1}_{n+1} )
        // + (omega^i_{n+1} - 1) . -------------------------------------------------
        //                                 | r^i_{n+1} - r^{i+1}_{n+1} |^2
        //
        // | r^i_{n+1} - r^{i+1}_{n+1} |^2
        // = | (-1)^2*(r^{i+1}_{n+1} - r^i_{n+1} |^2
        // = | (r^{i+1}_{n+1} - r^i_{n+1} |^2

        // initialise increment vector with solution of last iteration (i)
        // update del_ with current residual vector
        // difference of last two solutions
        if (del_ == Teuchos::null)
        {
         del_ = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
         delhist_ = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
         del_->PutScalar(1.0e20);
         delhist_->PutScalar(0.0);
        }

        // calculate difference of current (i+1) and old (i) residual vector
        // del = r^{i+1}_{n+1}
        del_->Update(1.0, *dispincnp_,0.0);
        // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
        delhist_->Update(1.0, *del_, (-1.0));
        double normdel = 0.0;
        double dot = 0.0;
        delhist_->Norm2(&normdel);
        // calculate dot product
        // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
        del_->Dot(*delhist_,&dot);

        if ((Step()==1) && (itnum == 1))
         dispnp->Update(1.0,*(StructureField()->Dispnp()),0.0);
        else  // (itnum > 1)
        {
         // relaxation parameter
         // omega^{i+1}_{n+1} == omeganp_
         // omega^{i}_{n+1} == omegan_
         // ome^{i+1} = ome^i + (ome^i -1) . (r^i - r^{i+1})^T . r^{i+1} / |r^{i+1} - r^{i}|^2
         omeganp_ = omegan_ + (omegan_ - 1.0) * (-dot) / (normdel*normdel);
         // relax displacement solution for next iteration step
         // overwrite dispnp with relaxed solution vector
         // d^{i+1} = d^{i+1} - omega^{i+1} * r^{i+1}
         //         = d^{i+1} - omega^{i+1} * ( d^{i+1} - d^i )
         //         = (1 - omega^{i+1}) d^{i+1} + omega^{i+1} d^i
         dispnp->Update(1.0,*(StructureField()->Dispnp()),(-omeganp_),*del_,0.0);
        }

        // update Aitken parameter omega^{i+1}_{n+1}
        omegan_ = omeganp_;

        // update history vector with residual displacement of old iteration step
        delhist_->Update(1.0,*del_,0.0);

        // end Aitken relaxation
        // ------------------------------------------------------------

        // update velocities
        // if relaxation active, then calculate new velocities using relaxed dispnp
        veln_ = CalcVelocity(dispnp);
        // TODO check if velocities have to be relaxed consistently to dispcacements

      } // end OUTER ITERATION

      // initial guess for next time step n+1
      // use minimum between omeganp_ and 1.0 as start value for omega_n+1^{i=0}
      // in case of doubt use 0.0, meaning that direction of new solution vector
      // dispnp better old one
      omegan_ = min(omeganp_, 0.0);

    }  // relax the displacements --> velocities
    else
      dserror("Relaxation of temperatures according to Irons & Tuck not available.");
  }  // end Aitken relaxation by Irons & Tuck

  return;
}  // OuterIterationLoop()


/*----------------------------------------------------------------------*
 | solve the structural system (protected)                    dano 12/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::DoStructureStep()
{
#ifndef TFSI
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }
#endif

  /// solve structural system
  // do the nonlinear solve for the time step. All boundary conditions have
  // been set.
  StructureField()->Solve();

  return;
}


/*----------------------------------------------------------------------*
 | solve the thermal system (protected)                       dano 12/10 |
 | coupling of displacements to thermo field and extract current        |
 | displacements needed for coupling to thermo field                    |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::DoThermoStep()
{
#ifndef TFSI
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }
#endif

  /// solve thermal system
  // do the solve for the time step. All boundary conditions have
  // been set.
  ThermoField()->Solve();

  return;
}


/*----------------------------------------------------------------------*
 | convergence check for both fields (thermo & structure)    dano 02/10 |
 | originally by vg 01/09                                               |
 *----------------------------------------------------------------------*/
bool TSI::Partitioned::ConvergenceCheck(
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

  // variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double tempincnorm_L2(0.0);
  double tempnorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1} with Newton iteration index i
  // \f Delta T^{i+1} = Inc T^{i+1} = T^{i+1} - T^{i}  \f
  tempincnp_->Update(1.0,*(ThermoField()->Tempnp()),-1.0);
  dispincnp_->Update(1.0,*(StructureField()->Dispnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  // for convergence test choose the last converged solution vector T_n/D_n,
  // be careful to check the convergence with the current, NOT yet converged values n+1
  // if a solution n+1 is highly difficult to find, the norm can oscillate
  ThermoField()->Tempn()->Norm2(&tempnorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  StructureField()->Dispn()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
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

  // converged
  // check convergence of the outer loop with a relative norm
  // norm of the increment with respect to the norm of the variables itself: use the last converged solution
  if ((tempincnorm_L2/tempnorm_L2 <= ittol) &&
      (dispincnorm_L2/dispnorm_L2 <= ittol))
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax);
      printf("+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum==itmax) and
       ((tempincnorm_L2/tempnorm_L2 > ittol) || (dispincnorm_L2/dispnorm_L2 > ittol))
     )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
    {
      printf("|     >>>>>> not converged in itemax steps!                                       |\n");
      printf("+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Partitioned::ConvergenceCheck


/*----------------------------------------------------------------------*/
