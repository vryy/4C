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
#include "tsi_defines.H"
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
#include "../drt_thermo/thr_contact.H"

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

  if (method==INPAR::TSI::OneWay)
  {
    if (displacementcoupling_) // (temperature change due to deformation)
    {
      // build a proxy of the structure discretization for the temperature field
      Teuchos::RCP<DRT::DofSet> structdofset
        = StructureField().Discretization()->GetDofSetProxy();
      // check if ThermoField has 2 discretizations, so that coupling is possible
      if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
        dserror("unexpected dof sets in thermo field");
    }

    else // temperature as coupling variable (deformation due to temperature change)
    {
      // build a proxy of the temperature discretization for the structure field
      Teuchos::RCP<DRT::DofSet> thermodofset
        = ThermoField().Discretization()->GetDofSetProxy();
      // check if StructField has 2 discretizations, so that coupling is possible
      if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
        dserror("unexpected dof sets in structure field");
    }
  }
  else if ( method==INPAR::TSI::SequStagg || method==INPAR::TSI::IterStagg )
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
  }

#ifdef TSIASOUTPUT
    // now check if the two dofmaps are available and then bye bye
    cout << "structure dofmap" << endl;
    cout << *StructureField().DofRowMap(0) << endl;
    cout << "thermo dofmap" << endl;
    cout << *StructureField().DofRowMap(1) << endl;
    cout << "thermo dofmap" << endl;
    cout << *ThermoField().DofRowMap(0) << endl;
    cout << "structure dofmap" << endl;
    cout << *ThermoField().DofRowMap(1) << endl;
//    exit(0);
#endif // TSIASOUTPUT

    // contact
    if(StructureField().ContactManager() != null)
      ThermoField().PrepareThermoContact(StructureField().ContactManager(),StructureField().Discretization());

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
 | update (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField().Update();
  ThermoField().Update();
  return;
}


/*----------------------------------------------------------------------*
 | output (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output()
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.
  StructureField().Output();

  ThermoField().Output();

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("UPRES");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");
  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if ( (upres!=0 and (Step()%upres == 0))
    or (uprestart != 0) and (Step()%uprestart == 0) )
    {
      // displacement field
      // (get noderowmap of discretisation for creating this multivector)
      dispnp_ = rcp(new Epetra_MultiVector(*(ThermoField().Discretization()->NodeRowMap()),3,true));

      OutputDeformationInThr(
          StructureField().Dispn(),
          StructureField().Discretization()
          );

      ThermoField().DiscWriter()->WriteVector("displacement",dispnp_,IO::DiscretizationWriter::nodevector);
    }

}

/*----------------------------------------------------------------------*
 | communicate the displacement vector to THR field          dano 12/11 |
 | enable visualisation of thermal variables on deformed body           |
 *----------------------------------------------------------------------*/
void  TSI::Algorithm::OutputDeformationInThr(
  Teuchos::RCP<const Epetra_Vector> dispnp,
  Teuchos::RCP<DRT::Discretization> structdis
  )
{
  if (dispnp == Teuchos::null)
    dserror("Got null pointer for displacements");

  int err(0);

  // get dofrowmap of structural discretisation
  const Epetra_Map* structdofrowmap = structdis->DofRowMap(0);

  // loop over all local nodes of thermal discretisation
  for (int lnodeid=0; lnodeid<(ThermoField().Discretization()->NumMyRowNodes()); lnodeid++)
  {
    // Here we rely on the fact that the thermal discretisation is a clone of
    // the structural mesh.
    // => a thermal node has the same local (and global) ID as its corresponding
    // structural node!

    // get the processor's local structural node with the same lnodeid
    DRT::Node* structlnode = structdis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this structural node
    std::vector<int> structnodedofs = structdis->Dof(0,structlnode);
    // determine number of space dimensions
    const int numdim = genprob.ndim;

    // now we transfer displacment dofs only
    for(int index=0; index<numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int sgid = structnodedofs[index];
      const int slid = structdofrowmap->LID(sgid);

      // get value of corresponding displacement component
      double disp = (*dispnp)[slid];
      // insert velocity value into node-based vector
      err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
      if (err!= 0) dserror("error while inserting a value into dispnp_");
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index=numdim; index < 3; ++index)
    {
      err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
      if (err!= 0) dserror("error while inserting a value into dispnp_");
    }

  } // for lnodid

return;

}  // OutputDeformationInThr()

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

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();

  // decide if apply one-way coupling or full coupling
  INPAR::TSI::SolutionSchemeOverFields method =
    DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // get active nodes from structural contact simulation
  RCP<MORTAR::ManagerBase> cmtman = StructureField().ContactManager();

  // tsi with or without contact
  // only tsi
  if (cmtman == Teuchos::null)
  {
    // time loop
    while (NotFinished())
    {
      // counter and print header
      IncrementTimeAndStep();
      PrintHeader();

      // normally PrepareTimeStep() should be called outside the nonlinear loop,
      // but at this point the coupling variables are not known, yet!
      // call PrepareTimeStep() after ApplyTemperatures()/ApplyStructVariables()

      // switch TSI algorithm (synchronous)
      // only the one field knows the 2nd field, that contains the coupling
      // variable
      if (method==INPAR::TSI::OneWay)
        TimeLoopOneWay();
      // complete volume coupling system
      // sequential staggered scheme
      else if (method==INPAR::TSI::SequStagg)
        TimeLoopSequStagg();
      // iterative staggered scheme
      else if (method==INPAR::TSI::IterStagg)
        TimeLoopFull();
      else
        dserror("desired type of thermo-structure interaction algorithm not supported");

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
      if(method==INPAR::TSI::OneWay)
      {
        // predict and solve structural system
        StructureField().PrepareTimeStep();
        StructureField().Solve();

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
      else if (method==INPAR::TSI::IterStagg)
      {
        // prepare time step
        ThermoField().PrepareTimeStep();
        StructureField().PrepareTimeStep();

        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;

        while (stopnonliniter==false)
        {
          itnum ++;

          // store temperature from first solution for convergence check (like in
          // elch_algorithm: use current values)
          tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
          dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);

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

      // calculate stresses, strains, energies
      PrepareOutput();

      // update all single field solvers
      Update();

      // write output to screen and files
      Output();
    } // time loop
  }

  // ==================================================================

}  // TSI::Algorithm::TimeLoop()


/*----------------------------------------------------------------------*
 | One-way coupling (only one field knows his counterpart)   dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopOneWay()
{
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    /// structure field

    // predict the structure field without influence of temperatures
    StructureField().PrepareTimeStep();

    /// solve structural system
    // get current displacement due to solve structure step
    DoStructureStep();

    // now extract the current displacements/velocities
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField().ExtractDispnp();

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
      velnp = StructureField().ExtractVelnp();
    }

    /// thermo field

    // use the structural solution calculated in DoStructureStep

    // call the current displacements and velocities and pass it to THR
    ThermoField().ApplyStructVariables(dispnp, velnp);

    // call the predictor here, because the displacement field is now also
    // available
    ThermoField().PrepareTimeStep();

    /// solve temperature system
    // do the solve for the time step. All boundary conditions have been set
    DoThermoStep();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // extract final displacements,
    // since we did update, this is very easy to extract
    dispn_ = StructureField().ExtractDispn();

  } // displacement coupling

  else // temperature coupling (deformation due to heating)
  {
    /// thermo field

    // predict the thermal field without influence of structure
    ThermoField().PrepareTimeStep();

    /// solve temperature system
    // get current temperatures due to Solve Thermo Step
    DoThermoStep();

    // now extract the current temperatures and pass it to the structure
    const Teuchos::RCP<Epetra_Vector> tempnp = ThermoField().ExtractTempnp();

    /// structure field

    // pass the current temperatures to the structural field
//    // 08.04.11
//    if( Step()==0 )
//      StructureField().ApplyTemperatures(tempnp,ThermoField().Tempn());
//    else
      StructureField().ApplyTemperatures(tempnp);

    // call the predictor here, because the temperature is now also available
    StructureField().PrepareTimeStep();

    /// solve structure system
    // do the nonlinear solve for the time step. All boundary conditions have
    // been set.
    DoStructureStep();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // extract final temperatures,
    // since we did update, this is very easy to extract
    tempn_ = ThermoField().ExtractTempn();

  }  // temperature coupling

}  // TSI::Algorithm::TimeLoopOneWay()


/*----------------------------------------------------------------------*
 | One-way coupling between the fields                       dano 09/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopSequStagg()
{
  // like implicit-implicit staggered scheme, compare Farhat & Park, 1991
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    // begin nonlinear solver ************************************************

    // predict the displacement field (structural predictor for TSI)

    // get structure variables of old time step (d_n, v_n)
    dispn_ = StructureField().ExtractDispn();
    if( Step()==0 ) veln_ = StructureField().ExtractVeln();

    /// thermo field

    // solve thermo field with predicted structural values (Apply d_n,v_n to
    // THR, solve coupled equation, extract new temperatures T_n+1)

    // pass the current displacements and velocities to the thermo field
    ThermoField().ApplyStructVariables(dispn_,veln_);

    // prepare time step with coupled variables
    ThermoField().PrepareTimeStep();

    /// solve temperature system

    // do the solve for the time step. All boundary conditions have been set.
    DoThermoStep();

    // now extract the current temperatures and pass it to the structure
    tempn_ = ThermoField().ExtractTempnp();

    /// structure field

    // pass the current temperatures to the structure field
    StructureField().ApplyTemperatures(tempn_);

    // prepare time step with coupled variables
    StructureField().PrepareTimeStep();

    // solve coupled equation
    DoStructureStep();

    // end nonlinear solver **************************************************

    // get the velocities needed as predictor for the thermo field for the next
    // time step
    if(quasistatic_)
      veln_ = CalcVelocity(StructureField().ExtractDispnp());
    else
      veln_ = StructureField().ExtractVelnp();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // extract final displacements,
    // since we did update, this is very easy to extract
    dispn_ = StructureField().ExtractDispn();

  } // end displacement coupling

  // temperature is coupling variable
  else  // (deformation due to heating)
  {
    // begin nonlinear solver ************************************************

    // predict the temperature field (thermal predictor for TSI)

    // get temperature solution of old time step
    tempn_ = ThermoField().ExtractTempn();

    /// structure field

    // get current displacement due to solve structure step

    // pass the current temperatures to the structure field
    StructureField().ApplyTemperatures(tempn_);

    // prepare time step with coupled variables
    StructureField().PrepareTimeStep();

    // solve structural coupled equation
    DoStructureStep();

    // extract current displacements
    const Teuchos::RCP<Epetra_Vector> dispnp = StructureField().ExtractDispnp();

    // extract the velocities of the current solution
    // the displacement -> velocity conversion
    if(quasistatic_)
      velnp_ = CalcVelocity(dispnp);
    // quasistatic to exlude oszillation
    else
      velnp_ = StructureField().ExtractVelnp();

    /// thermo field

    // pass the current displacements and velocities to the thermo field
    ThermoField().ApplyStructVariables(dispnp,velnp_);

    // prepare time step with coupled variables
    ThermoField().PrepareTimeStep();

    /// solve temperature system
    /// do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    DoThermoStep();

    // end nonlinear solver **************************************************

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // extract final displacements,
    // since we did update, this is very easy to extract
    tempn_ = ThermoField().ExtractTempn();

  } // temperature coupling

} // TSI::Algorithm::TimeLoopSequStagg()


/*----------------------------------------------------------------------*
 | Full coupling (iterative staggered scheme)                dano 06/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoopFull()
{
  // outer iteration loop
  OuterIterationLoop();

  // calculate stresses, strains, energies
  PrepareOutput();

  // update all single field solvers
  Update();

  // extract final displacement and velocity
  // since we did update, this is very easy to extract
  dispn_ = StructureField().ExtractDispn();
  tempn_ = ThermoField().ExtractTempn();

} // TSI::Algorithm::TimeLoopFull()


/*----------------------------------------------------------------------*
 | Outer Loop with convergence check                         dano 10/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::OuterIterationLoop()
{
  // iterate between the two fields
  int  itnum = 0;
  bool stopnonliniter = false;

  // outer iteration loop starts
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"**************************************************************\n";
    cout<<"      OUTER ITERATION LOOP \n";
    printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
      ThermoField().GetTimeNumStep());
    cout<<"**************************************************************\n";
  }

  // structural predictor
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    // mechanical predictor for the coupling iteration outside the loop
    // get structure variables of old time step (d_n, v_n)
    // d^p_n+1 = d_n, v^p_n+1 = v_n
    Teuchos::RCP<Epetra_Vector> dispnp;
    if ( Step()==1 )
    {
      dispnp = StructureField().ExtractDispn();
      veln_ = StructureField().ExtractVeln();
    }
    // else: use the velocity of the last converged step

    // start OUTER ITERATION
    while (stopnonliniter==false)
    {
      itnum ++;

      // store temperature from first solution for convergence check (like in
      // elch_algorithm: use current values)
      tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
      dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);

      // kind of mechanical predictor
      // 1st iteration: get structure variables of old time step (d_n, v_n)
      if (itnum==1) dispnp = StructureField().ExtractDispn();
      // for itnum>1 use the current solution dispnp of old iteration step

      // Begin Nonlinear Solver / Outer Iteration ******************************

      // thermo field

      // pass the current displacements and velocities to the thermo field
      ThermoField().ApplyStructVariables(dispnp,veln_);

      // prepare time step with coupled variables
      if (itnum==1)
        ThermoField().PrepareTimeStep();
      // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
      else if (itnum != 1)
        ThermoField().PreparePartitionStep();

      /// solve temperature system

      /// do the nonlinear solve for the time step. All boundary conditions have
      /// been set.
      DoThermoStep();

      // now extract the current temperatures and pass it to the structure
      tempn_ = ThermoField().ExtractTempnp();

      /// structure field

      // pass the current temperatures to the structure field
      // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
      // just calculated incremental solutions
      StructureField().ApplyTemperatures(tempn_);

      // prepare time step with coupled variables
      if (itnum==1)
        StructureField().PrepareTimeStep();
      // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
      else if (itnum != 1)
        StructureField().PreparePartitionStep();

      // solve coupled structural equation
      DoStructureStep();

      // extract current displacements
      const Teuchos::RCP<Epetra_Vector> dispnp = StructureField().ExtractDispnp();

      // end nonlinear solver / outer iteration ********************************

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

    // start OUTER ITERATION
    while (stopnonliniter==false)
    {
      itnum ++;
      // get current temperatures due to solve thermo step, like predictor in FSI
      if (itnum==1){ tempn_ = ThermoField().ExtractTempn(); }
      else { tempn_ = ThermoField().ExtractTempnp(); }

      // store temperature from first solution for convergence check (like in
      // elch_algorithm: here the current values are needed)
      tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);
      dispincnp_->Update(1.0,*StructureField().Dispnp(),0.0);

      // begin nonlinear solver / outer iteration ******************************

      /// structure field

      // pass the current temperatures to the structure field
      // ApplyTemperatures() also calls PreparePartitionStep() for prediction with
      // just calculated incremental solutions
      StructureField().ApplyTemperatures(tempn_);

      // prepare time step with coupled variables
      if (itnum==1)
        StructureField().PrepareTimeStep();
      // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
      else if (itnum != 1)
        StructureField().PreparePartitionStep();

      // solve coupled structural equation
      DoStructureStep();

      // Extract current displacements
      const Teuchos::RCP<Epetra_Vector> dispnp = StructureField().ExtractDispnp();

      // extract the velocities of the current solution
      // the displacement -> velocity conversion
      if(quasistatic_)
        velnp_ = CalcVelocity(dispnp);
      // quasistatic to exlude oszillation
      else
        velnp_ = StructureField().ExtractVelnp();

      /// thermo field

      // pass the current displacements and velocities to the thermo field
      ThermoField().ApplyStructVariables(dispnp,velnp_);

      // prepare time step with coupled variables
      if (itnum==1)
        ThermoField().PrepareTimeStep();
      // within the nonlinear loop, e.g. itnum>1 call only PreparePartitionStep
      else if (itnum != 1)
        ThermoField().PreparePartitionStep();

      /// solve coupled thermal system
      /// do the solve for the time step. All boundary conditions have been set.
      DoThermoStep();

      // end nonlinear solver / outer iteration ********************************

      // check convergence of both field for "partitioned scheme"
      stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

    } // end OUTER ITERATION

  } // temperature predictor

  return;
}  // OuterIterationLoop()


/*----------------------------------------------------------------------*
 | solve the structural system (protected)                    dano 12/10 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoStructureStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }

  /// solve structural system
  // do the nonlinear solve for the time step. All boundary conditions have
  // been set.
  StructureField().Solve();

  return;
}


/*----------------------------------------------------------------------*
 | solve the thermal system (protected)                       dano 12/10 |
 | coupling of displacements to thermo field and extract current        |
 | displacements needed for coupling to thermo field                    |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoThermoStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }

  /// solve thermal system
  // do the solve for the time step. All boundary conditions have
  // been set.
  ThermoField().Solve();

  return;
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

  // variables to save different L2 - Norms
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

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0)
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
  if ((tempincnorm_L2/tempnorm_L2 <= ittol) &&
      (dispincnorm_L2/dispnorm_L2 <= ittol))
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0)
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
    if ((Comm().MyPID()==0))
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
 | cf. InterfaceVelocity(disp) in FSI::DirichletNeumann                 |
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
}  // TSI::Algorithm::CalcVelocity


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
