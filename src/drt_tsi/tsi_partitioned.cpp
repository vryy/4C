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
#include "tsi_utils.H"
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
  tempincnp_(Teuchos::null),
  dispincnp_(Teuchos::null),
  disp_(Teuchos::null),
  veln_(Teuchos::null),
  velnp_(Teuchos::null),
  temp_(Teuchos::null),
  del_(Teuchos::null),
  delhist_(Teuchos::null),
  mu_(0.0)
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();

  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6
  normtypeinc_
    = DRT::INPUT::IntegralValue<INPAR::TSI::ConvNorm>(tsidyn,"NORM_INC");

  // decide if apply one-way coupling or full coupling
   coupling_ = DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // decide if one-way coupling or full coupling
//  INPAR::TSI::SolutionSchemeOverFields coupling =
//    DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");
  // coupling variable
  displacementcoupling_
    = tsidyn.get<std::string>("COUPVARIABLE") == "Displacement";
  if (displacementcoupling_)
    std::cout << "Coupling variable: displacement" << std::endl;
  else
    std::cout << "Coupling variable: temperature" << std::endl;

  // if structure field is quasi-static --> CalcVelocity
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // major switch to different time integrators
  quasistatic_
    = (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP") == INPAR::STR::dyna_statics);

  tempincnp_ = Teuchos::rcp(new Epetra_Vector(*(ThermoField()->Tempnp())));
  dispincnp_ = Teuchos::rcp(new Epetra_Vector(*(StructureField()->Dispnp())));
  temp_      = Teuchos::rcp(new Epetra_Vector(*(ThermoField()->Tempn())));

#ifdef TSIPARTITIONEDASOUTPUT
    // now check if the two dofmaps are available and then bye bye
    std::cout << "structure dofmap" << std::endl;
    std::cout << *StructureField()->DofRowMap(0) << std::endl;
    std::cout << "thermo dofmap" << std::endl;
    std::cout << *StructureField()->DofRowMap(1) << std::endl;
    std::cout << "thermo dofmap" << std::endl;
    std::cout << *ThermoField()->DofRowMap(0) << std::endl;
    std::cout << "structure dofmap" << std::endl;
    std::cout << *ThermoField()->DofRowMap(1) << std::endl;
//    exit(0);
#endif // TSIPARTITIONEDASOUTPUT

    // contact
    if (StructureField()->ContactManager() != Teuchos::null)
      ThermoField()->PrepareThermoContact(StructureField()->ContactManager(),StructureField()->Discretization());

}  // cstr


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

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  TSI::UTILS::SetMaterialPointersMatchingGrid(StructureField()->Discretization(),
                                              ThermoField()->Discretization());

  return;
}  // ReadRestart()


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
  disp_ = StructureField()->WriteAccessDispn();
  veln_ = StructureField()->WriteAccessVeln();

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | time loop                                                 dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::TimeLoop()
{

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
      // call PrepareTimeStep() after ApplyCouplingState()/ApplyStructVariables()

      // switch TSI algorithm (synchronous)
      // only the one field knows the 2nd field, that contains the coupling
      // variable
      // choose algorithm depending on solution type
      switch (coupling_)
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
      case INPAR::TSI::IterStaggFixedRel :
      {
        TimeLoopFull();
        break;
      }
      default :
        dserror("desired type of thermo-structure interaction algorithm not supported");
        break;
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
      if ((cmtman->GetStrategy().Friction()) == false)
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
      if (coupling_ == INPAR::TSI::OneWay)
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
      else if (coupling_ == INPAR::TSI::IterStagg)
      {
        // prepare time step
        ThermoField()->PrepareTimeStep();
        StructureField()->PrepareTimeStep();

        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;

        while (stopnonliniter == false)
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
          ApplyThermoCouplingState(tempnp);

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
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();

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
      velnp = StructureField()->WriteAccessVelnp();
    }

    /// thermo field

    // use the structural solution calculated in DoStructureStep

    // call the current displacements and velocities and pass it to THR
    ApplyStructCouplingState(dispnp, velnp);

    // call the predictor here, because the displacement field is now also
    // available
    ThermoField()->PrepareTimeStep();

    /// solve temperature system
    // do the solve for the time step. All boundary conditions have been set
    DoThermoStep();

    // extract final displacements,
    // since we did update, this is very easy to extract
    disp_ = StructureField()->WriteAccessDispnp();

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
    ApplyThermoCouplingState(tempnp);

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
    disp_ = StructureField()->WriteAccessDispn();
    if(Step() == 0) veln_ = StructureField()->WriteAccessVeln();

    /// thermo field

    // solve thermo field with predicted structural values (Apply d_n,v_n to
    // THR, solve coupled equation, extract new temperatures T_n+1)

    // pass the current displacements and velocities to the thermo field
    ApplyStructCouplingState(disp_,veln_);

    // prepare time step with coupled variables
    ThermoField()->PrepareTimeStep();

    /// solve temperature system

    // do the solve for the time step. All boundary conditions have been set.
    DoThermoStep();

    // now extract the current temperatures and pass it to the structure
    temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

    /// structure field

    // pass the current temperatures to the structure field
    ApplyThermoCouplingState(temp_);

    // prepare time step with coupled variables
    StructureField()->PrepareTimeStep();

    // solve coupled equation
    DoStructureStep();

    // end nonlinear solver **************************************************

    // get the velocities needed as predictor for the thermo field for the next
    // time step
    if(quasistatic_)
      veln_ = CalcVelocity(StructureField()->WriteAccessDispnp());
    else
      veln_ = StructureField()->WriteAccessVelnp();

    // extract final displacements,
    // update is called afterwards, so we extract the newest solution
    disp_ = StructureField()->WriteAccessDispnp();

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
    ApplyThermoCouplingState(temp_);

    // prepare time step with coupled variables
    StructureField()->PrepareTimeStep();

    // solve structural coupled equation
    DoStructureStep();

    // extract current displacements
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();

    // extract the velocities of the current solution
    // the displacement -> velocity conversion
    if(quasistatic_)
      velnp_ = CalcVelocity(dispnp);
    // quasistatic to exlude oszillation
    else
      velnp_ = StructureField()->WriteAccessVelnp();

    /// thermo field

    // pass the current displacements and velocities to the thermo field
    ApplyStructCouplingState(dispnp,velnp_);

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
  disp_ = StructureField()->WriteAccessDispnp();
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
  if ( (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
  {
    std::cout<<"\n";
    std::cout<<"**************************************************************\n";
    std::cout<<"      OUTER ITERATION LOOP \n";
    printf("      Time Step %3d/%3d \n",ThermoField()->GetTimeStep(),
      ThermoField()->GetTimeNumStep());
    std::cout<<"**************************************************************\n";
  }

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();
  // decide if one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(
        tsidyn,
        "COUPALGO"
        );

  // Pure iterative staggered algorithms
  // iterative staggered TSI withOUT Aitken relaxation
  if (coupling == INPAR::TSI::IterStagg)
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
      if (Step() == 1)
      {
        dispnp->Update(1.0, *(StructureField()->WriteAccessDispn()),0.0);
        veln_ = StructureField()->WriteAccessVeln();
      }
      // else: use the velocity of the last converged step

      // start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // kind of mechanical predictor
        // 1st iteration: get structure variables of old time step (d_n, v_n)
        if (itnum == 1)
          dispnp->Update(1.0,*(StructureField()->WriteAccessDispn()),0.0);
        // else (itnum>1) use the current solution dispnp of old iteration step

        // store temperature from first solution for convergence check (like in
        // elch_algorithm: use current values)
        tempincnp_->Update(1.0,*ThermoField()->Tempnp(),0.0);
        dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

        // begin nonlinear solver / outer iteration ****************************

        // ------------------------------------------------------- thermo field

        // pass the current displacements and velocities to the thermo field
        ApplyStructCouplingState(dispnp,veln_);

        // prepare time step with coupled variables
        if (itnum == 1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        // do the nonlinear solve for the time step. All boundary conditions
        // have been set.
        DoThermoStep();

        // now extract the current temperatures and pass it to the structure
        temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

        // ---------------------------------------------------- structure field

        // pass the current temperatures to the structure field
        //todo
        // ApplyCouplingState() also calls PreparePartitionStep() for prediction
        // with just calculated incremental solutions
        ApplyThermoCouplingState(temp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // extract current displacements
        dispnp->Update(1.0,*(StructureField()->WriteAccessDispnp()),0.0);

        if(quasistatic_)
          veln_ = CalcVelocity(dispnp);
        else
          veln_ = StructureField()->WriteAccessVelnp();

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        // end nonlinear solver / outer iteration ******************************

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
      while (stopnonliniter == false)
      {
        itnum ++;

        // get current temperatures due to solve thermo step, like predictor in FSI
        if (itnum == 1)
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

        // begin nonlinear solver / outer iteration ****************************

        // ---------------------------------------------------- structure field

        // pass the current temperatures to the structure field
        //todo
        // ApplyCouplingState() also calls PreparePartitionStep() for prediction
        // with just calculated incremental solutions
        ApplyThermoCouplingState(temp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // Extract current displacements
        const Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();

        // extract the velocities of the current solution
        // the displacement -> velocity conversion
        if(quasistatic_)
          velnp_ = CalcVelocity(dispnp);
        // quasistatic to exlude oszillation
        else
          velnp_ = StructureField()->WriteAccessVelnp();

        // ------------------------------------------------------- thermo field

        // pass the current displacements and velocities to the thermo field
        ApplyStructCouplingState(dispnp,velnp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        /// solve coupled thermal system
        /// do the solve for the time step. All boundary conditions have been set.
        DoThermoStep();

        // end nonlinear solver / outer iteration ******************************

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

      } // end OUTER ITERATION

    } // temperature predictor

  } // iterstagg WITHOUT relaxation

  // notation according to Aitken relaxation of FSI Mok, Uli
  // coupling==INPAR::TSI::IterStaggAitken
  // 1. calculate the Aitken factor mu
  // 2. calculate the relaxation factor omega = 1-mu
  // 3. T^{i+1} = T^i + omega^{i+1} * ( T^{i+1} - T^i )
  //            = T^{i+1} - mu_^{i+1} * ( T^{i+1} - T^i )
  // 4. limit Aitken factor mu_ for next time step with 1.0
  //
  // another notation of relaxation according to Paper by Irons & Tuck (1969)
  // coupling==INPAR::TSI::IterStaggAitkenIrons
  // 1. calculate an relaxation factor mu
  // 2. T^{i+1} = (1 - mu^{i+1}) T^{i+1} + mu^{i+1} T^i
  // 3. limit Aitken factor mu for next time step with maximal value MAXOMEGA
  //
  // coupling==INPAR::TSI::IterStaggFixedRel
  // 1. relaxation factor omega = FIXEDOMEGA = const
  // 2. T^{i+1} = omega^{i+1} . T^{i+1} + (1- omega^{i+1}) T^i
  //            = T^i + omega^{i+1} * ( T^{i+1} - T^i )
  else if ( (coupling == INPAR::TSI::IterStaggAitken) or
            (coupling == INPAR::TSI::IterStaggAitkenIrons) or
            (coupling == INPAR::TSI::IterStaggFixedRel)
          )
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << "Iterative staggered TSI with relaxation" << std::endl;
      std::cout << "Have you set MAXOMEGA to the wished value?" << std::endl;
    }

    // structural predictor
    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // mechanical predictor for the coupling iteration outside the loop
      // get structure variables of old time step (d_n, v_n)
      // d^p_n+1 = d_n, v^p_n+1 = v_n
      // initialise new time step n+1 with values of old time step n
      Teuchos::RCP<Epetra_Vector> dispnp
        = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
      if (Step() == 1)
      {
        dispnp->Update(1.0, *(StructureField()->WriteAccessDispn()),0.0);
        veln_ = StructureField()->WriteAccessVeln();
      }
      // else: use the velocity of the last converged step

      if ( (coupling == INPAR::TSI::IterStaggAitken) or
           (coupling == INPAR::TSI::IterStaggAitkenIrons)
         )
      {
        // constrain the Aitken factor in the 1st relaxation step of new time step
        // n+1 to maximal value maxomega
        double maxomega = tsidyn.get<double>("MAXOMEGA");
        // omega_{n+1} = min( omega_n, maxomega ) with omega = 1-mu
        if ( (maxomega > 0.0) and (maxomega < 1-mu_) )
          mu_ = 1 - maxomega;

        // reset in every new time step
        if (del_ != Teuchos::null)
        {
          del_->PutScalar(1.0e20);
        }
      }

      // start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // kind of mechanical predictor
        // 1st iteration: get structure variables of old time step (d_n, v_n)
        if (itnum == 1)
          dispnp->Update(1.0,*(StructureField()->WriteAccessDispn()),0.0);
        // else (itnum>1) use the current solution dispnp of old iteration step

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

        // begin nonlinear solver / outer iteration ****************************

        // ------------------------------------------------------- thermo field

        // pass the current displacements and velocities to the thermo field
        ApplyStructCouplingState(dispnp,veln_);

        // prepare time step with coupled variables
        if (itnum == 1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        // do the nonlinear solve for the time step. All boundary conditions
        // have been set.
        DoThermoStep();

        // now extract the current temperatures and pass it to the structure
        temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);

        // ---------------------------------------------------- structure field

        // pass the current temperatures to the structure field
        // ApplyCouplingState() also calls PreparePartitionStep() for prediction
        // with just calculated incremental solutions
        ApplyThermoCouplingState(temp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        if (coupling == INPAR::TSI::IterStaggFixedRel)
        {
          // get the displacements of the old iteration step d^i_{n+1}
          dispnp->Update(1.0,*(StructureField()->WriteAccessDispnp()),0.0);
        }

        // solve coupled structural equation
        DoStructureStep();

        // end nonlinear solver / outer iteration ******************************

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        // ------------------------------------------------------------
        // if r_{i+1} not converged in ConvergenceCheck()
        // --> apply relaxation to displacements

        if (coupling == INPAR::TSI::IterStaggFixedRel)
        {
          // -------------------------------------------- relax the displacements

          // get fixed relaxation parameter
          double fixedomega = tsidyn.get<double>("FIXEDOMEGA");
          // fixed relaxation can be applied even in the 1st iteration
          // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
          //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
          dispnp->Update(fixedomega,*dispincnp_,1.0);

          // -------------------------------------------------- end of relaxation
        }
        // dynamic relaxation
        else
        {
          // Aitken factor
          // relaxation can start at i=2, i.e. two vectors inc_{i-1}, inc_i are
          // available (also stated in Diss Uli, p. 82ff)
          // we need vectors from iteration step 1,2 to calculate relaxed step
          // with nu^i_0 = 0.0
          // Irons & Tuck:
          // increment r := old - new

          // BACI:
          // increment inc := new - old
          // nu^{i+1}_{n+1} = nu^i_{n+1} +
          //
          //                     ( r^{i+1}_{n+1} - r^i_{n+1} )^T . ( - r^{i+1}_{n+1} )
          // + (nu^i_{n+1} - 1) . -------------------------------------------------
          //                              | r^{i+1}_{n+1} - r^i_{n+1} |^2

          // initialise increment vector with solution of last iteration (i)
          // update del_ with current residual vector
          // difference of last two solutions
          if (del_ == Teuchos::null)  // first iteration, itnum==1
          {
            del_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
            delhist_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
            del_->PutScalar(1.0e20);
            delhist_->PutScalar(0.0);
          }

          // calculate difference of current (i+1) and old (i) residual vector
          // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
          // update history vector old increment r^i_{n+1}
          delhist_->Update(1.0,*del_,0.0);  // r^i_{n+1}
          delhist_->Update(1.0, *dispincnp_, (-1.0));  // update r^{i+1}_{n+1}

          // del_ = r^{i+1}_{n+1} = T^{i+1}_{n+1} - T^{i}_{n+1}
          del_->Update(1.0,*dispincnp_,0.0);
          // den = |r^{i+1} - r^{i}|^2
          double den = 0.0;
          delhist_->Norm2(&den);
          // calculate dot product
          // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
          double top = 0.0;
          delhist_->Dot(*del_,&top);

          // Aikten factor
          // nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2
          // top = ( r^{i+1} - r^i )^T . r^{i+1} --> use -top
          // Uli's implementation: mu_ = mu_ + (mu_ - 1.0) * top / (den*den). with '-' included in top
          mu_ = mu_ + (mu_ - 1)*(-top)/(den*den);

          // relaxation parameter
          // omega^{i+1} = 1- mu^{i+1}
          double omega = 1.0 - mu_;

          // two previous residuals are required (initial guess and result of 1st
          // iteration)
          // --> start relaxation process if two iterative residuals are available
          if (itnum == 1)
          {
            dispnp->Update(1.0,*(StructureField()->Dispnp()),0.0);

            // calculate the velocities with the updated/relaxed displacements
            if(quasistatic_)
              veln_ = CalcVelocity(dispnp);
            else
              veln_ = StructureField()->WriteAccessVelnp();
          }
          else  // (itnum > 1)
          {
            if (coupling == INPAR::TSI::IterStaggAitken)
            {
              // relax temperature solution for next iteration step
              // overwrite temp_ with relaxed solution vector
              // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
              //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
              dispnp->Update(omega,*dispincnp_,1.0);

            }
            // another notation of relaxation according to Paper by Irons & Tuck (1969)
            // 1. calculate an Aitken factor mu == relaxation factor
            // 2. d^{i+1} = (1 - mu^{i+1}) d^{i+1} + mu^{i+1} d^i
            // 3. limit Aitken factor mu for next time step with 0.0
            else if (coupling == INPAR::TSI::IterStaggAitkenIrons)
            {
              // relax displacement solution for next iteration step
              // overwrite dispnp with relaxed solution vector
              // d^{i+1} = d^{i+1} - mu^{i+1} * r^{i+1}
              //         = d^{i+1} - mu^{i+1} * ( d^{i+1} - d^i )
              //         = (1 - mu^{i+1}) d^{i+1} + mu^{i+1} d^i
              dispnp->Update(1.0,*(StructureField()->Dispnp()),(-mu_),*dispincnp_,0.0);
            }
          }  // itnum>1
        }  // dynamic relaxation

        // update velocities
        // velocity has to fit the corresponding displacements
        // if relaxation active --> calculate new velocities using relaxed dispnp
        veln_ = CalcVelocity(dispnp);

        // ----------------------------------------------------- end relaxation
      }  // end OUTER ITERATION
    }  // end relax displacement

    // relax the temperatures
    else // (!displacementcoupling_)
    {
      if (Step() == 1)
      {
        // thermal predictor for the coupling iteration outside the loop
        // get temperature of old time step (T_n)
        // T^p_n+1 = T_n
        temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
      }

      if ( (coupling == INPAR::TSI::IterStaggAitken) or
           (coupling == INPAR::TSI::IterStaggAitkenIrons)
         )
      {
        // initial guess for next time step n+1
        // use minimum between omega_n and maxomega as start value for mu_{n+1}^{i=0}
        // in case of doubt use 1.0, meaning that direction of new solution vector
        // tempnp better than old one

        // constrain the Aitken factor in the 1st relaxation step of new time step
        // n+1 to maximal value maxomega
        double maxomega = tsidyn.get<double>("MAXOMEGA");
        // omega_{n+1} = min( omega_n, maxomega )
        if ( (maxomega > 0.0) and (maxomega < 1-mu_) )
          mu_ = 1 - maxomega;

        // reset in every new time step
        if (del_ != Teuchos::null)
        {
          del_->PutScalar(1.0e20);
        }
      }

      // start OUTER ITERATION
      while (stopnonliniter == false)
      {
        itnum ++;

        // get current temperatures due to solve thermo step, like predictor in FSI
        // 1. iteration: get newest temperatures, i.e. T_{n+1} == T_n
        if (itnum == 1)
        {
          temp_->Update(1.0,*(ThermoField()->ExtractTempn()),0.0);
        }
        // else: use solution vector temp_ of old iteration step i

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

        // begin nonlinear solver / outer iteration ***************************

        // ---------------------------------------------------- structure field

        // pass the current temperatures to the structure field
        //todo
        // ApplyCouplingState() also calls PreparePartitionStep() for prediction
        // with just calculated incremental solutions
        //StructureField()->Discretization()->SetState(1,"temperature",temp_);
        ApplyThermoCouplingState(temp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          StructureField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          StructureField()->PreparePartitionStep();

        // solve coupled structural equation
        DoStructureStep();

        // Extract current displacements
        // now extract the current temperatures and pass it to the structure
        Teuchos::RCP<Epetra_Vector> dispnp
          = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
        dispnp->Update(1.0,*(StructureField()->WriteAccessDispnp()),0.0);

        // extract the velocities of the current solution
        // the displacement -> velocity conversion
        if(quasistatic_)
          velnp_ = CalcVelocity(dispnp);
        // quasistatic to exlude oscillation
        else
          velnp_ = StructureField()->WriteAccessVelnp();

        // ------------------------------------------------------- thermo field

        // pass the current displacements and velocities to the thermo field
        ApplyStructCouplingState(dispnp,velnp_);

        // prepare time step with coupled variables
        if (itnum == 1)
          ThermoField()->PrepareTimeStep();
        // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
        else if (itnum != 1)
          ThermoField()->PreparePartitionStep();

        if (coupling == INPAR::TSI::IterStaggFixedRel)
        {
          // get temperature solution of old iteration step T_{n+1}^i
          temp_->Update(1.0,*(ThermoField()->ExtractTempnp()),0.0);
        }

        /// solve coupled thermal system
        /// do the solve for the time step. All boundary conditions have been set.
        DoThermoStep();

        // end nonlinear solver / outer iteration ******************************

        // check convergence of both field for "partitioned scheme"
        stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);
        // --> update of tempincnp_ is done in ConvergenceCheck()

        // --------------------------------------------------- start relaxation

        if (coupling == INPAR::TSI::IterStaggFixedRel)
        {
          // ------------------------------------------- relax the temperatures

          // get fixed relaxation parameter
          double fixedomega = tsidyn.get<double>("FIXEDOMEGA");
          // fixed relaxation can be applied even in the 1st iteration
          // T^{i+1} = omega^{i+1} . T^{i+1} + (1- omega^{i+1}) T^i
          //         = T^i + omega^{i+1} * ( T^{i+1} - T^i )
          temp_->Update(fixedomega,*tempincnp_,1.0);

          // ------------------------------------------------ end of relaxation
        }
        // dynamic relaxation
        else
        {
          // ------------------------------------------------------------------
          // if r_{i+1} not converged in ConvergenceCheck()
          // --> apply Aitken relaxation to temperatures
          // as implemented in FSI by Irons and Tuck (1969)

          // Aitken factor
          // relaxation can start at i=2, i.e. two vectors inc_{i-1}, inc_i are
          // available (also stated in Diss Uli, p. 82ff)
          // we need vectors from iteration step 1,2 to calculate relaxed step
          // with mu^i_0 = 0.0
          // Irons & Tuck:
          // increment r := old - new

          // BACI:
          // increment inc := new - old
          // mu^{i+1}_{n+1} = mu^i_{n+1} +
          //
          //                     ( r^{i+1}_{n+1} - r^i_{n+1} )^T . ( - r^{i+1}_{n+1} )
          // + (mu^i_{n+1} - 1) . -------------------------------------------------
          //                              | r^{i+1}_{n+1} - r^i_{n+1} |^2

          // initialise increment vector with solution of last iteration (i)
          // update del_ with current residual vector
          // difference of last two solutions
          if (del_ == Teuchos::null)  // first iteration, itnum==1
          {
            del_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
            delhist_ = LINALG::CreateVector(*(ThermoField()->DofRowMap(0)), true);
            del_->PutScalar(1.0e20);
            delhist_->PutScalar(0.0);
          }

          // calculate difference of current (i+1) and old (i) residual vector
          // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
          // update history vector old increment r^i_{n+1}
          delhist_->Update(1.0,*del_,0.0);  // r^i_{n+1}
          delhist_->Update(1.0, *tempincnp_, (-1.0));  // update r^{i+1}_{n+1}

          // del_ = r^{i+1}_{n+1} = T^{i+1}_{n+1} - T^{i}_{n+1}
          del_->Update(1.0,*tempincnp_,0.0);
          // den = |r^{i+1} - r^{i}|^2
          double den = 0.0;
          delhist_->Norm2(&den);
          // calculate dot product
          // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
          double top = 0.0;
          delhist_->Dot(*del_,&top);

          // mu_: Aikten factor in Mok's version
          // mu_: relaxation parameter in Irons & Tuck
          // mu^{i+1} = mu^i + (mu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2
          // top = ( r^{i+1} - r^i )^T . r^{i+1} --> use -top
          // Uli's implementation: mu_ = mu_ + (mu_ - 1.0) * top / (den*den). with '-' included in top
          mu_ = mu_ + (mu_ - 1)*(-top)/(den*den);

          // two previous residuals are required (initial guess and result of 1st
          // iteration)
          // --> start relaxation process if two iterative residuals are available
          if (itnum == 1)
            // get temperature T_{n+1}^{i=1}
            temp_->Update(1.0,*(ThermoField()->Tempnp()),0.0);
          else  // (itnum > 1)
          {
            if (coupling == INPAR::TSI::IterStaggAitken)
            {
              // relaxation parameter
              // omega^{i+1} = 1- mu^{i+1}
              double omega = 1.0 - mu_;

              // relax temperature solution for next iteration step
              // overwrite temp_ with relaxed solution vector
              // T^{i+1} = omega^{i+1} . T^{i+1} + (1- omega^{i+1}) T^i
              //         = T^i + omega^{i+1} * ( T^{i+1} - T^i )
              temp_->Update(omega,*tempincnp_,1.0);
            }
            else if (coupling == INPAR::TSI::IterStaggAitkenIrons)
            {
              // relax displacement solution for next iteration step
              // overwrite tempnp with relaxed solution vector
              // T^{i+1} = T^{i+1} - mu^{i+1} * r^{i+1}
              //         = T^{i+1} - mu^{i+1} * ( T^{i+1} - T^i )
              //         = (1 - mu^{i+1}) T^{i+1} + mu^{i+1} T^i
              // --> Irons T^{i+1} = T^{i+1} + mu^{i+1} * DEL^{i+1} with DEL = T^i - T^{i+1}
              temp_->Update(1.0,*(ThermoField()->Tempnp()),(-mu_),*tempincnp_,0.0);
            }
          }  // itnum > 1
        }  // dynamic relaxation

        // ----------------------------------------------------- end relaxation

      } // end OUTER ITERATION
    } // relax temperatures
  } // iterative staggered TSI with relaxation
  return;
}  // OuterIterationLoop()


/*----------------------------------------------------------------------*
 | solve the structural system (protected)                    dano 12/10 |
 *----------------------------------------------------------------------*/
void TSI::Partitioned::DoStructureStep()
{
#ifndef TFSI
  if (Comm().MyPID() == 0)
  {
    std::cout<<"\n";
    std::cout<<"************************\n";
    std::cout<<"    STRUCTURE SOLVER    \n";
    std::cout<<"************************\n";
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
  if (Comm().MyPID() == 0)
  {
    std::cout<<"\n";
    std::cout<<"*********************\n";
    std::cout<<"    THERMO SOLVER    \n";
    std::cout<<"*********************\n";
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

  // converged
  // check convergence of the outer loop with a relative norm
  // norm of the increment with respect to the norm of the variables itself: use
  // the last converged solution
  // residual increments
  switch (normtypeinc_)
  {
    // default check:
    case INPAR::TSI::convnorm_abs:
    {
      // print the incremental based convergence check to the screen
      // test here increment
      if ( (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
      {
        std::cout<<"\n";
        std::cout<<"***********************************************************************************\n";
        std::cout<<"    OUTER ITERATION STEP    \n";
        std::cout<<"***********************************************************************************\n";
        printf("+--------------+------------------------+--------------------+--------------------+\n");
        printf("|-  step/max  -|-  tol      [norm]     -|--  temp-inc      --|--  disp-inc      --|\n");
        printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |",
        itnum,itmax,ittol,tempincnorm_L2,dispincnorm_L2);
        printf("\n");
        printf("+--------------+------------------------+--------------------+--------------------+\n");
      }

      // converged
      // check convergence of the outer loop with a relative norm
      // norm of the increment with respect to the norm of the variables itself: use the last converged solution
      stopnonliniter = ( (tempincnorm_L2 <= ittol) and (dispincnorm_L2 <= ittol) );
      if ( (stopnonliniter == true) and PrintScreenEvry() and
           (Comm().MyPID() == 0) and (Step()%PrintScreenEvry() == 0)
         )
      {
        printf("\n");
        printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax);
        printf("+--------------+------------------------+--------------------+--------------------+\n");
      }

      // warn if itemax is reached without convergence, but proceed to next
      // timestep
      if ( (itnum == itmax) and
           ( (tempincnorm_L2 > ittol) or (dispincnorm_L2 > ittol) )
         )
      {
        stopnonliniter = true;
        if ( (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
        {
          printf("|     >>>>>> not converged in itemax steps!                                       |\n");
          printf("+--------------+------------------------+--------------------+--------------------+\n");
          printf("\n");
          printf("\n");
        }
      }
    }  // INPAR::TSI::convnorm_abs
    break;

    case INPAR::TSI::convnorm_rel:
    {
      // print the incremental based convergence check to the screen
      // test here increment/variable
      if ( (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
      {
        std::cout<<"\n";
        std::cout<<"***********************************************************************************\n";
        std::cout<<"    OUTER ITERATION STEP    \n";
        std::cout<<"***********************************************************************************\n";
        printf("+--------------+------------------------+--------------------+--------------------+\n");
        printf("|-  step/max  -|-  tol      [norm]     -|--  temp-inc/temp --|--  disp-inc/disp --|\n");
        printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |",
        itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2,dispincnorm_L2/dispnorm_L2);
        printf("\n");
        printf("+--------------+------------------------+--------------------+--------------------+\n");
      }

      stopnonliniter = ( (tempincnorm_L2/tempnorm_L2 <= ittol) and (dispincnorm_L2/dispnorm_L2 <= ittol) );
      if ( (stopnonliniter == true) and (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
      {
        printf("\n");
        printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax);
        printf("+--------------+------------------------+--------------------+--------------------+\n");
      }

      // warn if itemax is reached without convergence, but proceed to next
      // timestep
      if ( (itnum == itmax) and
           ( (tempincnorm_L2/tempnorm_L2 > ittol) or (dispincnorm_L2/dispnorm_L2 > ittol) )
         )
      {
        stopnonliniter = true;
        if ( (Comm().MyPID() == 0) and PrintScreenEvry() and (Step()%PrintScreenEvry() == 0) )
        {
          printf("|     >>>>>> not converged in itemax steps!                                       |\n");
          printf("+--------------+------------------------+--------------------+--------------------+\n");
          printf("\n");
          printf("\n");
        }
      }
    }  // INPAR::TSI::convnorm_rel
    break;

    case INPAR::TSI::convnorm_mix:
    default:
      dserror("Cannot check for convergence of residual values!");
      break;
  }
  
  return stopnonliniter;
}  // TSI::Partitioned::ConvergenceCheck

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TSI::Partitioned::PrepareOutput()
{
  if(coupling_!=INPAR::TSI::OneWay or not displacementcoupling_)
    //set temperatures on structure field for evaluating stresses
    ApplyThermoCouplingState(ThermoField()->Tempnp());
  // prepare output (i.e. calculate stresses, strains, energies)
  StructureField()->PrepareOutput();

  //reset states
  StructureField()->Discretization()->ClearState(true);
}
/*----------------------------------------------------------------------*/
