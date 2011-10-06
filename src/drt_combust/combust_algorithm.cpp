/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

    detailed description in header file combust_algorithm.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_algorithm.H"
#include "combust_defines.H"
#include "two_phase_defines.H"
#include "combust_flamefront.H"
#include "combust_reinitializer.H"
#include "combust_utils.H"
#include "combust_fluidimplicitintegration.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_io/io_gmsh.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
  /* Der constructor sollte den gesamten Algorithmus initialisieren.
   * Hier müssen alle Variablen, die den Einzelfeldern übergeordnet sind, initialisiert werden.
   *
   * das heisst:
   * o G-Funktionsvektor (t und ig+1) auf initial value setzen
   * o Geschwindigkeitsvektor (t und iu+1) auf initial value setzen
   * o alle Zähler auf 0 setzen (step_(0), f_giter_(0), g_iter_(0), f_iter_(0))
   * o alle Normen und Grenzwerte auf 0 setzen
   *
   * Zusammenfassend kann an sagen, dass alles was in der Verbrennungsrechnung vor der Zeitschleife
   * passieren soll, hier passieren muss, weil die combust dyn gleich die Zeitschleife ruft.
   *
   * scalar transport velocity field has been initialized in ScaTraFluidCouplingAlgorithm()
  */
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn,false),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  convtol_(combustdyn.get<double>("CONVTOL")),
  stepbeforereinit_(false),
  stepreinit_(false),
  phireinitn_(Teuchos::null),
//  fgvelnormL2_(?),
//  fggfuncnormL2_(?),
  combusttype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  reinitaction_(DRT::INPUT::IntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
//  reinitaction_(combustdyn.sublist("COMBUSTION GFUNCTION").get<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION")),
  reinitinterval_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("REINITINTERVAL")),
  reinitband_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITBAND")),
  reinitbandwidth_(combustdyn.sublist("COMBUSTION GFUNCTION").get<double>("REINITBANDWIDTH")),
  reinit_output_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINIT_OUTPUT")),
  volcorrection_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITVOLCORRECTION")),
  extract_interface_vel_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"EXTRACT_INTERFACE_VEL")),
  convel_layers_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("NUM_CONVEL_LAYERS")),
  combustdyn_(combustdyn),
  interfacehandle_(Teuchos::null),
  interfacehandle_old_(Teuchos::null),
  flamefront_(Teuchos::null)
  {

  if (Comm().MyPID()==0)
  {
    switch(combusttype_)
    {
    case INPAR::COMBUST::combusttype_premixedcombustion:
      std::cout << "COMBUST::Algorithm: this is a premixed combustion problem" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem with kinks in vel and jumps in pres" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflowjump:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem with jumps in vel and pres" << std::endl;
      break;
    default: dserror("unknown type of combustion problem");
    }
  }

  // TODO:	scatra and fluid/combustion time-integration methods do not have to fit, or shall they? - winklmaier
  const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT");
  switch(timealgo)
  {
  case INPAR::FLUID::timeint_stationary:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
      dserror("fluid time integration scheme does not match");
    if (ScaTraField().MethodName() != INPAR::SCATRA::timeint_stationary)
      cout << "WARNING: combustion and scatra time integration scheme do not match" << endl;
    break;
  }
  case INPAR::FLUID::timeint_one_step_theta:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_one_step_theta)
      dserror("fluid time integration scheme does not match");
    if (ScaTraField().MethodName() != INPAR::SCATRA::timeint_one_step_theta)
      cout << "WARNING: combustion and scatra time integration scheme do not match" << endl;
    break;
  }
  case INPAR::FLUID::timeint_afgenalpha:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
      dserror("fluid time integration scheme does not match");
    if (ScaTraField().MethodName() != INPAR::SCATRA::timeint_gen_alpha)
      cout << "WARNING: combustion and scatra time integration scheme do not match" << endl;

    //TODO remove
    dserror("Generalized alpha time integration scheme for combustion is not working yet");
    break;
  }
  default:
    dserror("unknown time integration scheme for combustion problem");
  }

  // get pointers to the discretizations from the time integration scheme of each field
  // remark: fluiddis cannot be of type "const Teuchos::RCP<const DRT::Dis...>", because parent
  // class. InterfaceHandle only accepts "const Teuchos::RCP<DRT::Dis...>"              henke 01/09
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  velnpip_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);
  velnpi_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);

  // FGI vectors are initialized
  velnpip_->Update(1.0,*FluidField().Velnp(),0.0);
  velnpi_->Update(1.0,*FluidField().Velnp(),0.0);

  phinpip_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
  phinpi_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);

  // FGI vectors are initialized
  phinpip_->Update(1.0,*ScaTraField().Phinp(),0.0);
  phinpi_->Update(1.0,*ScaTraField().Phinp(),0.0);

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   *
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/
  // construct initial flame front
  flamefront_ = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis,ScaTraField().PBCmap()));
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  // construct interfacehandles using initial flame front
  interfacehandle_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));
  interfacehandle_old_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));
  // get integration cells according to initial flame front
  interfacehandle_->UpdateInterfaceHandle();
  interfacehandle_old_->UpdateInterfaceHandle();

  volume_start_ = ComputeVolume(); // schott

  if (reinitaction_ != INPAR::COMBUST::reinitaction_none)
  {
//    stepreinit_ = true;
    DoReinitialization(); // schott 04/11

//	  cout<< "reinitialization at timestep 0 switched off!" << endl;

    if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary)
    {
      // reset phin vector in ScaTra time integration scheme to phinp vector
      *ScaTraField().Phin() = *ScaTraField().Phinp();
    }
    // pointer not needed any more
    stepreinit_ = false;
  }

  //------------------------
  // set initial fluid field
  //------------------------
  const INPAR::COMBUST::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(
      combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD");
  const int initfuncno = combustdyn_.sublist("COMBUSTION FLUID").get<int>("INITFUNCNO");

  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);

  FluidField().SetInitialFlowField(initfield, initfuncno);

  // output fluid initial state
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary)
    FluidField().Output();

  // TODO check if needed
  // export interface information to the fluid time integration
  // remark: this is essential here, if DoFluidField() is not called in Timeloop() (e.g. for pure Scatra problems)
  FluidField().ImportInterface(interfacehandle_,interfacehandle_old_);

  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion)
  {
    // extract convection velocity from fluid solution
    const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();
    ScaTraField().SetVelocityField(
//        OverwriteFluidVel(),
        //FluidField().ExtractInterfaceVeln(),
        ComputeFlameVel(convel,FluidField().DofSet()),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
    );
  }

  // output G-function initial state
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary and
      DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false )
    ScaTraField().Output();
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a dynamic combustion problem                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{
  // compute initial volume of minus domain
  volume_start_ = ComputeVolume();

  // get initial field by solving stationary problem first
  if(DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == true)
    SolveInitialStationaryProblem();

  // time loop
  while (NotFinished())
  {
    // prepare next time step; update field vectors
    PrepareTimeStep();

    // Fluid-G-function-Interaction loop
    while (NotConvergedFGI())
    {
      // prepare Fluid-G-function iteration
      PrepareFGIteration();

      // TODO In the first iteration of the first time step the convection velocity for the
      //      G-function is zero, if a zero initial fluid field is used.
      //      -> Should the fluid be solved first?
      // solve linear G-function equation
      DoGfuncField();

      //(after Scatra transport but before reinitialization)
      // update interface geometry
      UpdateInterface();

      // TODO: stepreinit and fg-iteration (Benedikt)
      // if reinitialization is accepted the interface and interfacehandle was updated
      if(stepreinit_) DoReinitialization();


      // solve nonlinear Navier-Stokes system
      DoFluidField();

    } // Fluid-G-function-Interaction loop

    // update all field solvers
    UpdateTimeStep();
    //Remark (important for restart): the time level of phi (n+1, n or n-1) used to reconstruct the interface
    //                                conforming to the restart state of the fluid depends on the order
    //                                of Output() and UpdateTimeStep()

    // write output to screen and files
    Output();

    if (!stepreinit_)
    {
      // compute current volume of minus domain
      const double volume_current = ComputeVolume();
      // print mass conservation check on screen
      printMassConservationCheck(volume_start_, volume_current);
    }

  } // time loop

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start_, volume_end);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary combustion problem                              henke 10/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("--------Stationary-Combustion-------  time step ----------------------------------------\n");
  }
  //-----------------------------
  // prepare stationary algorithm
  //-----------------------------
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // check if ScaTraField().initialvelset == true
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField().MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // compute initial volume of minus domain
  const double volume_start = ComputeVolume();

  //--------------------------------------
  // loop over fluid and G-function fields
  //--------------------------------------
  while (NotConvergedFGI())
  {
    // prepare Fluid-G-function iteration
    PrepareFGIteration();

    // solve nonlinear Navier-Stokes system
    DoFluidField();

    // solve (non)linear G-function equation
    std::cout << "/!\\ warning === G-function field not solved for stationary problems" << std::endl;
    //DoGfuncField();

    //reinitialize G-function
    //if (fgiter_ % reinitinterval_ == 0)
    //{
    //  ReinitializeGfunc();
    //}

    // update field vectors
    UpdateInterface();

  } // fluid-G-function loop

  //-------
  // output
  //-------
  // remark: if Output() was already called at initial state, another Output() call will cause an
  //         error, because both times fields are written into the output control file at time and
  //         time step 0.
  //      -> the time and the time step have to be advanced even though this makes no physical sense
  //         for a stationary computation
  //IncrementTimeAndStep();
  //FluidField().PrepareTimeStep();
  //ScaTraField().PrepareTimeStep();

  // write output to screen and files (and Gmsh)
  Output();

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start, volume_end);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                schott 04/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareReinitialization()
{

	// reset the ScatraFieldReinit
	ScaTraReinitField().ResetTimeAndStep();

	// set the start phinp and phin field and phistart field
	ScaTraReinitField().SetPhiReinit(ScaTraField().Phinp());

	// set params
	// set velocity field for reinitialization -> implement this in Scatra-part



	return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                schott 05/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoReinitialization()
{
//  if(reinitaction_ == INPAR::COMBUST::reinitaction_pdebased_characteristic_galerkin)
//	  {cout << "Reinitialization Characteristic "
//	 reinitaction_ == INPAR::COMBUST::reinitaction_pdebased_stabilized_convection)


  // compute current volume of minus domain
  const double volume_current_before = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start_, volume_current_before);


  // reinitialize Gfunc
  switch(reinitaction_)
  {
    case INPAR::COMBUST::reinitaction_sussman:
      if(Comm().MyPID()==0) std::cout << "COMBUST::Algorithm: reinitialization via PDE-based method" << std::endl;
      PrepareReinitialization();
      reinitialization_accepted_ = ScaTraReinitField().CallReinitialization();
      if(reinitialization_accepted_ == true) ScaTraField().SetPhinp(ScaTraReinitField().Phinp());
      break;
    case INPAR::COMBUST::reinitaction_signeddistancefunction:
    case INPAR::COMBUST::reinitaction_fastsigneddistancefunction:
      if(Comm().MyPID()==0) std::cout << "COMBUST::Algorithm: reinitialization via calculating signed distance functions" << std::endl;
      ReinitializeGfuncSignedDistance();
      reinitialization_accepted_ = true;
      break;
    case INPAR::COMBUST::reinitaction_none:
      if(Comm().MyPID()==0) std::cout << "No reinitialization chosen" << std::endl;
      reinitialization_accepted_=false;
      break;
    default: dserror("unknown type of reinitialization technique");
  }

  // now ScatraField().Phinp() is updated

  if(Comm().MyPID()==0)
  {
    if(reinitialization_accepted_) cout << "...reinitialization was accepted." << std::endl;
    else                           cout << "...reinitialization was not accepted." << std::endl;
  }



  // update flamefront and interfacehandle according to reinitialized G-function field
  if(reinitialization_accepted_)
  {
    // update flame front according to reinitialized G-function field
    flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

    // update interfacehandle (get integration cells) according to updated flame front
    interfacehandle_->UpdateInterfaceHandle();
  }

  // TODO: the step used for output is switched for one
  if(reinit_output_ and reinitialization_accepted_ and !(reinitaction_ == INPAR::COMBUST::reinitaction_signeddistancefunction or reinitaction_ == INPAR::COMBUST::reinitaction_fastsigneddistancefunction))
  {
    ScaTraReinitField().OutputReinit(ScaTraField().Step(), ScaTraField().Time());
  }

  // compute current volume of minus domain
  const double volume_current_after = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start_, volume_current_after);

  // do volume correction
  if (volcorrection_)
  {
    CorrectVolume(volume_start_, volume_current_after);

    // compute current volume of minus domain
    const double volume_current_corrected = ComputeVolume();
    // print mass conservation check on screen
    printMassConservationCheck(volume_start_, volume_current_corrected);
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                schott 04/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfuncSignedDistance()
{
	// REMARK:
	// all interface updates and updates of interfacehandles are done in DoReinitialization() (schott)

    // get my flame front (boundary integration cells)
    // TODO we generate a copy of the flame front here, which is not neccessary in the serial case
    std::map<int, GEO::BoundaryIntCells> myflamefront = interfacehandle_->GetElementalBoundaryIntCells();

#ifdef PARALLEL
    // export flame front (boundary integration cells) to all processors
    flamefront_->ExportFlameFront(myflamefront);
#endif

    // reinitialize G-function (level set) field
    COMBUST::Reinitializer reinitializer(
        combustdyn_,
        ScaTraField(),
        myflamefront,
        ScaTraField().Phinp());

    return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc()
{
  //  // turn on/off screen output for writing process of Gmsh postprocessing file
  //  const bool screen_out = true;
  //  // create Gmsh postprocessing file
  //  const std::string filename1 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_reinitialization_before", Step(), 701, screen_out, ScaTraField().Discretization()->Comm().MyPID());
  //  std::ofstream gmshfilecontent1(filename1.c_str());
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent1 << "View \" " << "Phinp \" {" << endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(ScaTraField().Discretization(),ScaTraField().Phinp(),gmshfilecontent1);
  //    gmshfilecontent1 << "};" << endl;
  //  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent1 << "View \" " << "Convective Velocity \" {" << endl;
  //    // draw vector field 'Convective Velocity' for every element
  //    IO::GMSH::VectorFieldNodeBasedToGmsh(ScaTraField().Discretization(),ScaTraField().ConVel(),gmshfilecontent1);
  //    gmshfilecontent1 << "};" << endl;
  //  }
  //  gmshfilecontent1.close();
  //  if (screen_out) std::cout << " done" << endl;

  // reinitialize what will later be 'phinp'
  if (stepreinit_)
  {
    // REMARK:
    // maybe there are some modified phi-values in phinp_ and so the current interfacehandle_ is adopted
    // to the modified values
    // we want to reinitialize the orignial phi-values
    // => UpdateFlameFront without modifyPhiVectors

    // update flame front according to reinitialized G-function field
    // the reinitilizer needs the original G-function field
    // ModifyPhiVector uses an alternative modification such that tetgen does not crash, so UpdateFlameFront is called
    // with the boolian true

    if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") == INPAR::FLUID::timeint_afgenalpha)
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phiaf());
    else
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());

    // update interfacehandle (get integration cells) according to updated flame front
    interfacehandle_->UpdateInterfaceHandle();


    // get my flame front (boundary integration cells)
    // TODO we generate a copy of the flame front here, which is not neccessary in the serial case
    std::map<int, GEO::BoundaryIntCells> myflamefront = interfacehandle_->GetElementalBoundaryIntCells();

#ifdef PARALLEL
    // export flame front (boundary integration cells) to all processors
    flamefront_->ExportFlameFront(myflamefront);
#endif

    if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") == INPAR::FLUID::timeint_afgenalpha)
    {
      // reinitialize G-function (level set) field
      COMBUST::Reinitializer reinitializer(
          combustdyn_,
          ScaTraField(),
          myflamefront,
          ScaTraField().Phiaf());

      // after the reinitialization we update the flamefront in the usual sense
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phiaf());
    }
    else
    {
      // reinitialize G-function (level set) field
      COMBUST::Reinitializer reinitializer(
          combustdyn_,
          ScaTraField(),
          myflamefront,
          ScaTraField().Phinp());

      // after the reinitialization we update the flamefront in the usual sense
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());
    }

    // update interfacehandle (get integration cells) according to updated flame front
    interfacehandle_->UpdateInterfaceHandle();
  }

  //  // create Gmsh postprocessing file
  //  const std::string filename2 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_reinitialization_after", Step(), 701, screen_out, ScaTraField().Discretization()->Comm().MyPID());
  //  std::ofstream gmshfilecontent2(filename2.c_str());
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent2 << "View \" " << "Phinp \" {" << endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(ScaTraField().Discretization(),ScaTraField().Phinp(),gmshfilecontent2);
  //    gmshfilecontent2 << "};" << endl;
  //  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent2 << "View \" " << "Convective Velocity \" {" << endl;
  //    // draw vector field 'Convective Velocity' for every element
  //    IO::GMSH::VectorFieldNodeBasedToGmsh(ScaTraField().Discretization(),ScaTraField().ConVel(),gmshfilecontent2);
  //    gmshfilecontent2 << "};" << endl;
  //  }
  //  gmshfilecontent2.close();
  //  if (screen_out) std::cout << " done" << endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: overwrite Navier-Stokes velocity                                        henke 08/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::OverwriteFluidVel()
{
  //----------------------------------------------------------------
  // impose velocity field by function e.g. for level set test cases
  // Navier-Stokes solution velocity field is overwritten
  //----------------------------------------------------------------
  if (Comm().MyPID()==0)
    std::cout << "\n--- overwriting Navier-Stokes solution with field defined by FUNCT1... " << std::flush;

  // get fluid (Navier-Stokes) velocity vector in standard FEM configuration (no XFEM dofs)
  const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();
  // velocity function number = 1 (FUNCT1)
  const int velfuncno = 1;

  // loop all nodes on the processor
  for(int lnodeid=0; lnodeid < FluidField().Discretization()->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    DRT::Node*  lnode = FluidField().Discretization()->lRowNode(lnodeid);
    // get standard dofset from fluid time integration
    vector<int> fluidnodedofs = (*(FluidField().DofSet())).Dof(lnode);
    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;
    if (numdim != 3) dserror("3 components expected for velocity");

    int err = 0;
    // overwrite velocity dofs only
    for(int icomp=0; icomp<numdim; ++icomp)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[icomp];
      int flid = convel->Map().LID(fgid);
      if (flid < 0) dserror("lid not found in map for given gid");

      // get value of corresponding velocity component
      double value = DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(icomp,lnode->X(),0.0,NULL);

//      // scaling for stretching fluid example // schott
//      cout <<"Scaling with time-curve!!!!" << endl;
//      value *= cos(PI*FluidField().Time()/500.0);
      //value = 0.0;

      // insert velocity value into node-based vector
      err += convel->ReplaceMyValues(1, &value, &flid);
    }
    if (err != 0) dserror("error overwriting Navier-Stokes solution");
  }

  if (Comm().MyPID()==0)
    std::cout << "done" << std::endl;

  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | protected: compute reinit velocity                                                schott 04/11 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ComputeReinitVel(
    /*const Teuchos::RCP<Epetra_Vector>& convel,*/
    const Teuchos::RCP<const DRT::DofSet>& dofset
    )
{
	// needed for lids
	const Teuchos::RCP<Epetra_Vector> & convel = FluidField().ExtractInterfaceVeln();

	  // temporary vector for convective velocity (based on dofrowmap of standard (non-XFEM) dofset)
	  // remark: operations must not be performed on 'convel', because the vector is accessed by both
	  //         master and slave nodes, if periodic bounday conditions are present
	  Teuchos::RCP<Epetra_Vector> reinitvel_tmp = LINALG::CreateVector(*(dofset->DofRowMap()),true);

	  // get a pointer to the fluid discretization
	  const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();

	  // get G-function value vector on fluid NodeColMap
	  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
	  const Teuchos::RCP<Epetra_MultiVector> gradphinp = flamefront_->GradPhi();

	#ifdef DEBUG
	  // get map of this vector
	  const Epetra_BlockMap& phimap = phinp->Map();
	  // check, whether this map is still identical with the current node map in the discretization
	  if (not phimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
	  // get map of this vector
	  const Epetra_BlockMap& gradphimap = gradphinp->Map();
	  // check, whether this map is still identical with the current node map in the discretization
	  if (not gradphimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
	#endif

	#ifdef COMBUST_GMSH_NORMALFIELD
	  const std::string filestr = "normal_field";
	  const std::string name_in_gmsh = "Normal field";

	  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 500, true, fluiddis->Comm().MyPID());
	  std::ofstream gmshfilecontent(filename.c_str());
	  {
	    gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
	#endif

	#if 1
	    // loop over nodes on this processor
	    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
	    {
	      //-------------------------------------------------------------
	      // get (smoothed) gradient of the G-function field at this node
	      //-------------------------------------------------------------
	      LINALG::Matrix<3,1> mygradphi(true);

	      // get the processor local node
	      DRT::Node*  lnode = fluiddis->lRowNode(lnodeid);
	      // get the global id for current node
	      const int gid = lnode->Id();
	      // get local processor id according to global node id
	      const int nodelid = (*gradphinp).Map().LID(gid);
	      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphinp).Comm().MyPID(),gid);

	      const int numcol = (*gradphinp).NumVectors();
	      if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not 3");

	      // loop over dimensions (= number of columns in multivector)
	      for(int col=0; col< numcol; col++)
	      {
	        // get columns vector of multivector
	        double* globalcolumn = (*gradphinp)[col];
	        // set smoothed gradient entry of phi into column of global multivector
	        mygradphi(col) = globalcolumn[nodelid];
	      }

	      //------------------------------------
	      // compute smoothed unit normal vector
	      // n = - grad phi / |grad phi|
	      //------------------------------------
	      // smoothed normal vector at this node
	      LINALG::Matrix<3,1> nvec = mygradphi;
	      // compute norm of smoothed normal vector
	      const double ngradphi = mygradphi.Norm2();

	      // divide vector by its norm to get unit normal vector
	      if (abs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
	      {
	        // length of smoothed normal is zero at this node -> node must be on a singularity of the
	        // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
	        // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
	        //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
	        std::cout << "/!\\ reinit velocity at node " << gid << " is set to zero" << std::endl;
	        nvec.Clear();
	      }
	      else
	      {
	    	  nvec.Scale(1.0/ngradphi);
	      }

	#ifdef COMBUST_GMSH_NORMALFIELD
	      LINALG::SerialDenseMatrix xyz(3,1);
	      xyz(0,0) = lnode->X()[0];
	      xyz(1,0) = lnode->X()[1];
	      xyz(2,0) = lnode->X()[2];

	      IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, nvec, xyz, gmshfilecontent);
	#endif

	      //------------------------
//	      // get material parameters
//	      //------------------------
//	      // get list of adjacent elements of this node
//	      DRT::Element** elelist = lnode->Elements();
//
//	      // get material from first (arbitrary!) element adjacent to this node
//	      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
//	      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
//	      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());
//
//	      // density burnt domain
//	      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
//	      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
//	      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
//	      const double rhoplus = matplus->Density();
//
//	      // density unburnt domain
//	      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
//	      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
//	      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
//	      const double rhominus = matminus->Density();
//
//	      // laminar flame speed
//	      const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
//	      //---------------------------------------------
//	      // compute relative flame velocity at this node
//	      //---------------------------------------------
//	      // get phi value for this node
//	      const double gfuncval = (*phinp)[nodelid];
//
//	      double speedfac = 0.0;
//	      if (gfuncval >= 0.0){ // interface or burnt domain -> burnt material
//	        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
//	        speedfac = sl * rhominus/rhoplus;
//	      }
//	      else{ // unburnt domain -> unburnt material
//	        // flame speed factor = laminar flame speed
//	        speedfac = sl;
//	      }


//	      cout << "WARNING: the velocity vector is set to zero at the moment" << endl;
//	      nvec.Clear();

	      //-----------------------------------------------
	      // compute (absolute) flame velocity at this node
	      //-----------------------------------------------
	      LINALG::Matrix<3,1> fluidvel(true);
	      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
	      const std::vector<int> dofids = (*dofset).Dof(lnode);
	      std::vector<int> lids(3);

	      // extract velocity values (no pressure!) from global velocity vector
	      for (int icomp=0; icomp<3; ++icomp)
	      {
	        lids[icomp] = convel->Map().LID(dofids[icomp]);
	      }

	      LINALG::Matrix<3,1> flvelabs(true);
	      // add fluid velocity (Navier Stokes solution) and relative flame velocity
	      for (int icomp=0; icomp<3; ++icomp)
	      {
	        const int err = reinitvel_tmp->ReplaceMyValues(1,&nvec(icomp),&lids[icomp]);
	        if (err) dserror("could not replace values for convective reinit velocity");
	      }
	    }
	#endif

	#ifdef COMBUST_GMSH_NORMALFIELD
	    gmshfilecontent << "};\n";
	  }
	  gmshfilecontent.close();
	  if (true) std::cout << " done" << endl;
	#endif

  return reinitvel_tmp;
}

/*------------------------------------------------------------------------------------------------*
 | protected: compute flame velocity                                                  henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ComputeFlameVel(
    const Teuchos::RCP<Epetra_Vector>& convel,
    const Teuchos::RCP<const DRT::DofSet>& dofset
    )
{
  if((DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false) and
      (DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD") == INPAR::COMBUST::initfield_zero_field))
    cout << "/!\\ Compute an initial stationary fluid solution to avoid a non-zero initial flame velocity" << endl;

  // temporary vector for convective velocity (based on dofrowmap of standard (non-XFEM) dofset)
  // remark: operations must not be performed on 'convel', because the vector is accessed by both
  //         master and slave nodes, if periodic bounday conditions are present
  Teuchos::RCP<Epetra_Vector> conveltmp = LINALG::CreateVector(*(dofset->DofRowMap()),true);

  // get a pointer to the fluid discretization
  const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();

  // get G-function value vector on fluid NodeColMap
  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
  const Teuchos::RCP<Epetra_MultiVector> gradphinp = flamefront_->GradPhi();
  const Teuchos::RCP<Epetra_Vector> curvature = flamefront_->Curvature();

#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
  // get map of this vector
  const Epetra_BlockMap& gradphimap = gradphinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not gradphimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
#endif

  // laminar flame speed
  const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
  const double D = 0.0;//1.161/(1.161*1.0);

#ifdef COMBUST_GMSH_NORMALFIELD
  const std::string filestr = "normal_field";
  const std::string name_in_gmsh = "Normal field";

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 500, true, fluiddis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
#endif

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      //-------------------------------------------------------------
      // get (smoothed) gradient of the G-function field at this node
      //-------------------------------------------------------------
      LINALG::Matrix<3,1> mygradphi(true);

      // get the processor local node
      DRT::Node*  lnode = fluiddis->lRowNode(lnodeid);
      // get the global id for current node
      const int gid = lnode->Id();
      // get local processor id according to global node id
      const int nodelid = (*gradphinp).Map().LID(gid);
      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphinp).Comm().MyPID(),gid);

      const int numcol = (*gradphinp).NumVectors();
      if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not 3");

      // loop over dimensions (= number of columns in multivector)
      for(int col=0; col< numcol; col++)
      {
        // get columns vector of multivector
        double* globalcolumn = (*gradphinp)[col];
        // set smoothed gradient entry of phi into column of global multivector
        mygradphi(col) = globalcolumn[nodelid];
      }

      //------------------------------------
      // compute smoothed unit normal vector
      // n = - grad phi / |grad phi|
      //------------------------------------
      // smoothed normal vector at this node
      LINALG::Matrix<3,1> nvec = mygradphi;
      // compute norm of smoothed normal vector
      const double ngradphi = mygradphi.Norm2();

      // divide vector by its norm to get unit normal vector
      if (abs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
      {
        // length of smoothed normal is zero at this node -> node must be on a singularity of the
        // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
        // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
        //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
        std::cout << "/!\\ phi gradient too small at node " << gid << " -> flame velocity is only the convective velocity" << std::endl;
        nvec.PutScalar(0.0);
      }
      else
      {
        nvec.Scale(-1.0/ngradphi);
      }

#ifdef COMBUST_GMSH_NORMALFIELD
      LINALG::SerialDenseMatrix xyz(3,1);
      xyz(0,0) = lnode->X()[0];
      xyz(1,0) = lnode->X()[1];
      xyz(2,0) = lnode->X()[2];

      IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, nvec, xyz, gmshfilecontent);
#endif

      //------------------------
      // get material parameters
      //------------------------
      // get list of adjacent elements of this node
      DRT::Element** elelist = lnode->Elements();

      // get material from first (arbitrary!) element adjacent to this node
      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

      // density burnt domain
      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
      const double rhoplus = matplus->Density();

      // density unburnt domain
      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
      const double rhominus = matminus->Density();

      //---------------------------------------------
      // compute relative flame velocity at this node
      //---------------------------------------------
      // get phi value for this node
      const double gfuncval = (*phinp)[nodelid];
      const double curv = 0.0;//(*curvature)[nodelid];

      double speedfac = 0.0;
      if (gfuncval >= 0.0){ // interface or burnt domain -> burnt material
        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
        speedfac = (sl -D*curv)* rhominus/rhoplus;
      }
      else{ // unburnt domain -> unburnt material
        // flame speed factor = laminar flame speed
        speedfac = sl-D*curv;
      }

      LINALG::Matrix<3,1> flvelrel(true);
      for (int icomp=0; icomp<3; ++icomp)
        flvelrel(icomp) = speedfac * nvec(icomp);

      //-----------------------------------------------
      // compute (absolute) flame velocity at this node
      //-----------------------------------------------
      LINALG::Matrix<3,1> fluidvel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);
      std::vector<int> lids(3);

      // extract velocity values (no pressure!) from global velocity vector
      for (int icomp=0; icomp<3; ++icomp)
      {
        lids[icomp] = convel->Map().LID(dofids[icomp]);
        fluidvel(icomp) = (*convel)[lids[icomp]];
      }

      LINALG::Matrix<3,1> flvelabs(true);
      // add fluid velocity (Navier Stokes solution) and relative flame velocity
      for (int icomp=0; icomp<3; ++icomp)
      {
        flvelabs(icomp) = fluidvel(icomp) + flvelrel(icomp);
        const int err = conveltmp->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
        if (err) dserror("could not replace values for convective velocity");
      }
    }
#ifdef COMBUST_GMSH_NORMALFIELD
    gmshfilecontent << "};\n";
  }
  // Gmsh output curvature
  {
    gmshfilecontent << "View \" " << "Curvature \" {\n";
    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode = fluiddis->lRowNode(lnodeid);
      const int gid = lnode->Id();
      const int nodelid = (*curvature).Map().LID(gid);
      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphinp).Comm().MyPID(),gid);
      const double curv = (*curvature)[nodelid];
      LINALG::SerialDenseMatrix xyz(3,1);
      xyz(0,0) = lnode->X()[0];
      xyz(1,0) = lnode->X()[1];
      xyz(2,0) = lnode->X()[2];

      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, curv, xyz, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  gmshfilecontent.close();
  if (true) std::cout << " done" << endl;
#endif

  return conveltmp;
}

/*------------------------------------------------------------------------------------------------*
 | protected: FGI iteration converged?                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
  // if (fgiter <= fgitermax and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
  //return (fgiter_ < fgitermax_ and true);

  bool notconverged = true;
  
  // TODO just a first ordinary nonparallel error estimator
  // implement this parallel and for velocity if fgi-iter > 1 is used

  // TODO remove as soon as FG iteration works
    phinpi_->Update(1.0,*(ScaTraField().Phinp()),0.0);
#if 0
/* at the moment only the convergence of the g-function field is checked
 * to check the convergence of the fluid field, uncomment the corresponding lines
 */
  if (fgiter_ == 0)
  {
    // G-function field and fluid field haven't been solved, yet

    // store old solution vectors
//    if (velnpip_->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
//      dserror("vectors must have the same length 1");
//    velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
    phinpi_->Update(1.0,*(ScaTraField().Phinp()),0.0);

    if (fgiter_ == fgitermax_)
    {
      notconverged = false;
      if (Comm().MyPID()==0)
        cout << "!!!WARNING: A maximum of 0 FGI is not sensible!!!" << endl;

    }
  }
  else if (fgiter_ > 0)
  {
//    double velnormL2 = 1.0;
    double gfuncnormL2 = 1.0;

    // store new solution vectors and compute L2-norm
//    velnpi_->Update(1.0,*velnpip_,0.0);
    phinpi_->Update(1.0,*phinpip_,0.0);
//    velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
//    velnpip_->Norm2(&velnormL2);
    phinpip_->Update(1.0,*(ScaTraField().Phinp()),0.0);
    phinpip_->Norm2(&gfuncnormL2);

//    if (velnormL2 < 1e-5) velnormL2 = 1.0;
    if (gfuncnormL2 < 1e-5) gfuncnormL2 = 1.0;

//    fgvelnormL2_ = 0.0;
    fggfuncnormL2_ = 0.0; //TODO warum diese norm als klassenattribut und nicht lokal nur hier?

    // compute increment and L2-norm of increment
//    Teuchos::RCP<Epetra_Vector> incvel = rcp(new Epetra_Vector(velnpip_->Map()),true);
//    if (incvel->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
//      dserror("vectors must have the same length 2");
//    incvel->Update(1.0,*velnpip_,-1.0,*velnpi_,0.0);
//    incvel->Norm2(&fgvelnormL2_);
    Teuchos::RCP<Epetra_Vector> incgfunc = rcp(new Epetra_Vector(*ScaTraField().Discretization()->DofRowMap()),true);
    incgfunc->Update(1.0,*phinpip_,-1.0,*phinpi_,0.0);
    incgfunc->Norm2(&fggfuncnormL2_);

    if (fgitermax_ > 1)
    // and if (Comm().MyPID()==0) // TODO ein, wenn parallel geht!
    {
      cout << "on proc " << Comm().MyPID() << endl;
      printf("\n|+------------------------ FGI ------------------------+|");
      printf("\n|iter/itermax|----tol-[Norm]--|-fluid inc--|-g-func inc-|");
      printf("\n|   %2d/%2d    | %10.3E[L2] | ---------- | %10.3E |",fgiter_,fgitermax_,convtol_,fggfuncnormL2_/gfuncnormL2);
      printf("\n|+-----------------------------------------------------+|\n");

      convtol_ = 1e-16; // TODO remove if works...
      if (fggfuncnormL2_/gfuncnormL2 <= convtol_) //((fgvelnormL2_/velnormL2 <= convtol_) and (fggfuncnormL2_/gfuncnormL2 <= convtol_))
      {
        notconverged = false;
      }
      else
      {
        if (fgiter_ == fgitermax_)
        {
          notconverged = false;
          if (Comm().MyPID()==0)
          {
            printf("|+---------------- not converged ----------------------+|");
            printf("\n|+-----------------------------------------------------+|\n");
          }
        }
      }
    }
  }
#endif

  if (fgiter_ >= fgitermax_)
    notconverged = false;

  return notconverged;
}

/*------------------------------------------------------------------------------------------------*
 | protected: do a stationary first time step for combustion algorithm               schott 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveInitialStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("==============================================================================================\n");
    printf("----------------Stationary timestep prepares instationary algorithm---------------------------\n");
    printf("==============================================================================================\n");
  }
  //-----------------------------
  // prepare stationary algorithm
  //-----------------------------
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // check if ScaTraField().initialvelset == true
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // modify time and timestep for stationary timestep
  SetTimeStep(0.0,0); // algorithm timestep

  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }

  FluidField().PrepareTimeStep();
  
  // compute initial volume of minus domain
  const double volume_start = ComputeVolume();

  //-------------------------------------
  // solve nonlinear Navier-Stokes system
  //-------------------------------------
  DoFluidField();

  // update field vectors
  UpdateInterface();

  // assign the fluid velocity field to the G-function as convective velocity field
  switch(combusttype_)
  {
  case INPAR::COMBUST::combusttype_twophaseflow:
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    // for two-phase flow, the fluid velocity field is continuous; it can be directly transferred to
    // the scalar transport field

    ScaTraField().SetVelocityField(
      //OverwriteFluidVel(),
      FluidField().ExtractInterfaceVeln(),
      Teuchos::null,
      FluidField().DofSet(),
      FluidField().Discretization()
    );

    // Transfer history vector only for subgrid-velocity
    //ScaTraField().SetVelocityField(
    //    FluidField().ExtractInterfaceVeln(),
    //    FluidField().Hist(),
    //    FluidField().DofSet(),
    //    FluidField().Discretization()
    //);
    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added

    // extract convection velocity from fluid solution
    const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();

    ScaTraField().SetVelocityField(
//        OverwriteFluidVel(),
        //FluidField().ExtractInterfaceVeln(),
        ComputeFlameVel(convel,FluidField().DofSet()),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
    );
    break;
  }
  default:
    dserror("unknown type of combustion problem");
  }


  //-------
  // output
  //-------
  // remark: if Output() was already called at initial state, another Output() call will cause an
  //         error, because both times fields are written into the output control file at time and
  //         time step 0.
  //      -> the time and the time step have to be advanced even though this makes no physical sense
  //         for a stationary computation
  //IncrementTimeAndStep();
  //FluidField().PrepareTimeStep();
  //ScaTraField().PrepareTimeStep();

  // write output to screen and files (and Gmsh)
  Output();

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start, volume_end);
  
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // TODO @Florian clarify if this parameter is still needed
//  stepbeforereinit_ = false;
//  if (Step()>0 and Step() % reinitinterval_ == 0) stepbeforereinit_ = true;
  stepreinit_ = false;
  if (Step() % reinitinterval_ == 0) stepreinit_ = true; //Step()>1 and

  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }

  FluidField().PrepareTimeStep();

  // prepare time step
  // remark: initial velocity field has been transferred to scalar transport field in constructor of
  //         ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such
  //         as the one-step-theta scheme, are thus initialized correctly.
  ScaTraField().PrepareTimeStep();

  // synchronicity check between combust algorithm and base algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid time integration differs from time in combustion algorithm");
  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra time integration  differs from time in combustion algorithm");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
  fgiter_ += 1;
  if (Comm().MyPID()==0)
  {
    //cout<<"\n---------------------------------------  FGI loop  -------------------------------------------\n";
    printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
  }
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  if (Comm().MyPID()==0)
  {
    std:: cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------" << std::endl;
  }

  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  // export interface information to the fluid time integration
  FluidField().ImportInterface(interfacehandle_,interfacehandle_old_);
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  // solve nonlinear Navier-Stokes equations
  FluidField().NonlinearSolve();
  //FluidField().MultiCorrector();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a G-function time integration step                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncField()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  G-FUNCTION SOLVER  ----------------------------------\n";
  }
  // assign the fluid velocity field to the G-function as convective velocity field
  switch(combusttype_)
  {
  case INPAR::COMBUST::combusttype_twophaseflow:
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    // for two-phase flow, the fluid velocity field is continuous; it can be directly transferred to
    // the scalar transport field

    if (extract_interface_vel_) // replace computed velocity field by the interface velocity
    {
      ScaTraField().SetVelocityField(
        //OverwriteFluidVel(),
        ManipulateFluidFieldForGfunc(FluidField().ExtractInterfaceVeln(), FluidField().DofSet()),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
      );
    }
    else // transfer the computed velocity as convective velocity to the g-function solver
    {
      ScaTraField().SetVelocityField(
        //OverwriteFluidVel(),
        FluidField().ExtractInterfaceVeln(),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
      );

      // Transfer history vector only for subgrid-velocity
      //ScaTraField().SetVelocityField(
      //    FluidField().ExtractInterfaceVeln(),
      //    FluidField().Hist(),
      //    FluidField().DofSet(),
      //    FluidField().Discretization()
      //);
    }

    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added

    // extract convection velocity from fluid solution
    const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();

#if 0
    //std::cout << "convective velocity is transferred to ScaTra" << std::endl;
    Epetra_Vector convelcopy = *convel;
    //    *copyconvel = *convel;

    const Teuchos::RCP<Epetra_Vector> xfemvel = ComputeFlameVel(convel,FluidField().DofSet());
    if((convelcopy).MyLength() != (*xfemvel).MyLength())
      dserror("length is not the same!");

    const int dim = (convelcopy).MyLength();
    int counter = 0;

    for(int idof=0; idof < dim ;++idof)
    {
      if ((convelcopy)[idof] == (*xfemvel)[idof])
        counter++;
    }
    // number of identical dofs in convection velocity and flame velocity vector
    cout << "number of identical velocity components: " << counter << endl;
#endif

    ScaTraField().SetVelocityField(
//        OverwriteFluidVel(),
        //FluidField().ExtractInterfaceVeln(),
        ComputeFlameVel(convel,FluidField().DofSet()),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
    );
    break;
  }
  default:
    dserror("unknown type of combustion problem");
  }

  //solve convection-diffusion equation
  ScaTraField().Solve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateInterface()
{
  //overwrite old interfacehandle before updating flamefront
  interfacehandle_old_->UpdateInterfaceHandle();

  // update flame front according to evolved G-function field
  // remark: for only one FGI iteration, 'phinpip_' == ScaTraField().Phin()
  // TODO @Martin Bitte checken wie das hier korrekt geht, so ist es bestimmt falsch!
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") == INPAR::FLUID::timeint_afgenalpha)
    flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phiaf());
  else
    flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());
//  flamefront_->UpdateFlameFront(combustdyn_,phinpip_, ScaTraField().Phinp()); // phinpip_ now is the secondary newest phi-vector

  // update interfacehandle (get integration cells) according to updated flame front
  interfacehandle_->UpdateInterfaceHandle();

  // update the Fluid and the FGI vector at the end of the FGI loop
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField().Update();

  if (stepreinit_ and reinitialization_accepted_)
  {
	cout << "UpdateReinit" << endl;
    ScaTraField().UpdateReinit();
  }
  else
  {
     cout << "Update" << endl;
     ScaTraField().Update();
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  //------------------------------------------------------------------------------------------------
  // this hack is necessary for the visualization of disconituities in Gmsh             henke 10/09
  //------------------------------------------------------------------------------------------------
  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  FluidField().Output();
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  // causes error in DEBUG mode (trueresidual_ is null)
  //FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


void COMBUST::Algorithm::printMassConservationCheck(const double volume_start, const double volume_end)
{
  if (Comm().MyPID() == 0)
  {
    // compute mass loss
    if (volume_start != 0.0)
    {

      double massloss = -(volume_start - volume_end) / volume_start *100;
      // 'isnan' seems to work not reliably; error occurs in line above
      if (std::isnan(massloss))
        dserror("NaN detected in mass conservation check");

      std::cout << "---------------------------------------" << endl;
      std::cout << "           mass conservation           " << endl;
      std::cout << " initial mass: " << volume_start << endl;
      std::cout << " final mass:   " << volume_end   << endl;
      std::cout << " mass loss:    " << massloss << "%" << endl;
      std::cout << "---------------------------------------" << endl;
    }
    else
    {
      dserror(" there is no 'minus domain'! -> division by zero checking mass conservation");
    }
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | compute volume on all processors                                                   henke 02/10 |
 *------------------------------------------------------------------------------------------------*/
double COMBUST::Algorithm::ComputeVolume()
{
  // compute negative volume of discretization on this processor
  double myvolume = interfacehandle_->ComputeVolumeMinus();

  double sumvolume = 0.0;

  // sum volumes on all processors
  // remark: ifndef PARALLEL sumvolume = myvolume
  Comm().SumAll(&myvolume,&sumvolume,1);

  return sumvolume;
}


/*------------------------------------------------------------------------------------------------*
 | correct the volume of the minus domain after reinitialization                  rasthofer 07/11 |
 |                                                                                    DA wichmann |
 | Idea: shift level-set so that volume is conserved                                              |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::CorrectVolume(const double targetvol, const double currentvol)
{
  const double voldelta = targetvol - currentvol;

  if (Comm().MyPID() == 0)
    cout << "Correcting volume of minus(-) domain by " << voldelta << " ... " << flush;

  // compute negative volume of discretization on this processor
  double myarea = interfacehandle_->ComputeSurface();
  double sumarea = 0.0;

  cout << "Proc " << Comm().MyPID() << " myarea " << myarea << endl;

  // sum volumes on all processors
  // remark: ifndef PARALLEL sumvolume = myvolume
  Comm().SumAll(&myarea,&sumarea,1);

  if (Comm().MyPID() == 0)
    cout << "sumarea " << sumarea << endl;

  // This is a guess on how thick a layer needs to be added to the surface of the minus domain.
  // Due to $ \grad \phi \approx 1 $ this also happens to be the value that needs to be subtracted
  // of all phis. To make sure that \grad \phi really is close to 1 this function should only be
  // called after a reinitialization.
  const double thickness = -voldelta / sumarea;

  if (Comm().MyPID() == 0)
    cout << "thickness " << thickness << endl;

  RCP<Epetra_Vector> phin = ScaTraField().Phin();

  RCP<Epetra_Vector> one = rcp(new Epetra_Vector(phin->Map()));
  one->PutScalar(1.0);

  if (phin != Teuchos::null)
    phin->Update(thickness, *one, 1.0);

  RCP<Epetra_Vector> phinp = ScaTraField().Phinp();
  if (phinp != Teuchos::null)
    phinp->Update(thickness, *one, 1.0);

  RCP<Epetra_Vector> phinm = ScaTraField().Phinm();
  if (phinm != Teuchos::null)
    phinm->Update(thickness, *one, 1.0);

  // REMARK:
  // after the reinitialization we update the flamefront in the usual sense
  // this means that we modify phi-values if necessary -> default boolian modifyPhiVectors = true

  // update flame front according to reinitialized G-function field
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  // update interfacehandle (get integration cells) according to updated flame front
  interfacehandle_->UpdateInterfaceHandle();

  if (Comm().MyPID() == 0)
    cout << "done" << endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 |                                                                                rasthofer 08/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ManipulateFluidFieldForGfunc(
    const Teuchos::RCP<Epetra_Vector>& convel,
    const Teuchos::RCP<const DRT::DofSet>& dofset
    )
{
  if (convel_layers_ > 0)
  {
    // temporary vector for convective velocity (based on dofrowmap of standard (non-XFEM) dofset)
    // remark: operations must not be performed on 'convel', because the vector is accessed by both
    //         master and slave nodes, if periodic bounday conditions are present
    Teuchos::RCP<Epetra_Vector> conveltmp = LINALG::CreateVector(*(dofset->DofRowMap()),true);

    // get a pointer to the fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();

    // get G-function value vector on fluid NodeColMap
    const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();

    // some variables necessary for the Gather command
    //const int thisproc = fluiddis->Comm().MyPID();
    const int numproc  = fluiddis->Comm().NumProc();
    int allproc[numproc];
    for (int i = 0; i < numproc; ++i)
      allproc[i] = i;

    //--------------------------------------------------------------------------------------------------
    // due to the PBCs we need to get some info here in order to properly handle it later
    //--------------------------------------------------------------------------------------------------
    Teuchos::RCP<PeriodicBoundaryConditions> pbc = ScaTraField().PBC();

    // get the following information about the pbc
    // - planenormaldirection e.g. (1,0,0)
    // - minimum in planenormaldirection
    // - maximum in planenormaldirection
    vector<DRT::Condition*>* surfacepbcs = pbc->ReturnSurfacePBCs();
    vector<int>    planenormal(0);
    vector<double> globalmins (0);
    vector<double> globalmaxs (0);
    for (size_t i = 0; i < surfacepbcs->size(); ++i)
    {
      const string* ismaster = (*surfacepbcs)[i]->Get<string>("Is slave periodic boundary condition");
      if (*ismaster == "Master")
      {
        const int masterid = (*surfacepbcs)[i]->GetInt("Id of periodic boundary condition");
        vector<int> nodeids(*((*surfacepbcs)[i]->Nodes()));
        for (size_t j = 0; j < surfacepbcs->size(); ++j)
        {
          const int slaveid = (*surfacepbcs)[j]->GetInt("Id of periodic boundary condition");
          if (masterid == slaveid)
          {
            const string* isslave = (*surfacepbcs)[j]->Get<string>("Is slave periodic boundary condition");
            if (*isslave == "Slave")
            {
              const vector<int>* slavenodeids = (*surfacepbcs)[j]->Nodes();
              // append slave node Ids to node Ids for the complete condition
              for (size_t k = 0; k < slavenodeids->size(); ++k)
               nodeids.push_back(slavenodeids->at(k));
            }
          }
        }

        // Get normal direction of pbc plane
        const string* pbcplane = (*surfacepbcs)[i]->Get<string>("degrees of freedom for the pbc plane");
        if (*pbcplane == "yz")
          planenormal.push_back(0);
        else if (*pbcplane == "xz")
          planenormal.push_back(1);
        else if (*pbcplane == "xy")
          planenormal.push_back(2);
        else
          dserror("A PBC condition could not provide a plane normal.");

        double min = +10e19;
        double max = -10e19;
        for (size_t j = 0; j < nodeids.size(); ++j)
        {

          const int gid = nodeids[j];
          const int lid = fluiddis->NodeRowMap()->LID(gid);
          if (lid < 0)
            continue;
          const DRT::Node* lnode = fluiddis->lRowNode(lid);
          const double* coord = lnode->X();
          if (coord[planenormal.back()] < min)
            min = coord[planenormal.back()];
          if (coord[planenormal.back()] > max)
           max = coord[planenormal.back()];
        }
        globalmins.resize(planenormal.size());
       globalmaxs.resize(planenormal.size());
        fluiddis->Comm().MinAll(&min, &(globalmins.back()), 1);
        fluiddis->Comm().MaxAll(&max, &(globalmaxs.back()), 1);
      }
    }//end loop over all surfacepbcs


    // these sets contain the element/node GIDs that have been collected
    Teuchos::RCP<set<int> > allcollectednodes    = rcp(new set<int>);
    Teuchos::RCP<set<int> > allcollectedelements = rcp(new set<int>);

    // this loop determines how many layers around the cut elements will be collected
    for (int loopcounter = 0; loopcounter < convel_layers_; ++loopcounter)
    {
      if (loopcounter == 0)
      {
        //-------------------------------------------------------------------------------------------------
        // loop over row elements an check whether it has positive and negative phi-values. If it does
        // add the element to the allcollectedelements set.
        //-------------------------------------------------------------------------------------------------
        for(int lroweleid=0; lroweleid < fluiddis->NumMyRowElements(); ++lroweleid)
        {
          DRT::Element* ele = fluiddis->lRowElement(lroweleid);
          const int* nodeids = ele->NodeIds();
          bool gotpositivephi = false;
          bool gotnegativephi = false;

          for(int inode = 0; inode < ele->NumNode(); ++inode)
          {
            const int nodegid = nodeids[inode];
            const int nodelid = phinp->Map().LID(nodegid);
            if (nodelid < 0)
              dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*phinp).Comm().MyPID(),nodegid);

            if ((*phinp)[nodelid] <= 0.0)
              gotnegativephi = true;
            else
              gotpositivephi = true;
          }

          if (gotpositivephi and gotnegativephi)
            allcollectedelements->insert(ele->Id());
        }
      }
      else
      {
        //------------------------------------------------------------------------------------------
        // now, that we have collected all row nodes for this proc. Its time to get the adjacent elements
        //------------------------------------------------------------------------------------------
        set<int>::const_iterator nodeit;
        for (nodeit = allcollectednodes->begin(); nodeit != allcollectednodes->end(); ++nodeit)
        {
          const int nodelid       = fluiddis->NodeRowMap()->LID(*nodeit);
          DRT::Node* node         = fluiddis->lRowNode(nodelid);
          DRT::Element** elements = node->Elements();
          for (int iele = 0; iele < node->NumElement(); ++iele)
          {
            DRT::Element* ele     = elements[iele];
            allcollectedelements->insert(ele->Id());
          }
        } // loop over elements
      }

      //------------------------------------------------------------------------------------------------
      // now that we have collected all elements on this proc, its time to get the adjacent nodes.
      // Afterwards the nodes are communicated in order to obtain the collected row nodes on every proc.
      //------------------------------------------------------------------------------------------------
      set<int>::const_iterator eleit;
      for (eleit = allcollectedelements->begin(); eleit != allcollectedelements->end(); ++eleit)
      {
        const int elelid  = fluiddis->ElementColMap()->LID(*eleit);
        DRT::Element* ele = fluiddis->lColElement(elelid);
        DRT::Node** nodes = ele->Nodes();
        for (int inode = 0; inode < ele->NumNode(); ++inode)
        {
          DRT::Node* node = nodes[inode];

          // now check whether we have a pbc condition on this node
          vector<DRT::Condition*> mypbc;
          node->GetCondition("SurfacePeriodic",mypbc);

          if (mypbc.size() == 0)
          {
            allcollectednodes->insert(node->Id());
          }
          else
          {
            // obtain a vector of master and slaves
            const int nodeid = node->Id();
            vector<int> pbcnodes;
            for (size_t numcond=0; numcond<mypbc.size(); ++numcond)
            {
              Teuchos::RCP<map<int,vector<int> > > pbcmastertoslave = ScaTraField().PBCmap();
              map<int,vector<int> >::iterator mapit;
              for (mapit = pbcmastertoslave->begin(); mapit != pbcmastertoslave->end(); ++mapit)
              {
                if (mapit->first == nodeid)
                {
                  pbcnodes.push_back(mapit->first);
                  for (size_t i = 0; i < mapit->second.size(); ++i)
                    pbcnodes.push_back(mapit->second[i]);
                  break;
                }
                for (size_t isec = 0; isec < mapit->second.size(); ++isec)
                {
                  if (mapit->second[isec] == nodeid)
                  {
                    pbcnodes.push_back(mapit->first);
                    for (size_t i = 0; i < mapit->second.size(); ++i)
                      pbcnodes.push_back(mapit->second[i]);
                    break;
                  }
                }
              }
            }

            for (size_t i = 0; i < pbcnodes.size(); ++i)
              allcollectednodes->insert(pbcnodes[i]);
          }
        } // loop over elements' nodes
      } // loop over elements

      // with all nodes collected it is time to communicate them to all other procs
      // which then eliminate all but their row nodes
      {
        Teuchos::RCP<set<int> > globalcollectednodes = rcp(new set<int>);
        //TODO Ursula: patch gfunctionConvelManipulation
        //             linalg_utils.H, drt_exporter.H
        LINALG::Gather<int>(*allcollectednodes,*globalcollectednodes,numproc,allproc,fluiddis->Comm());

        allcollectednodes->clear();
        set<int>::const_iterator gnodesit;
        for (gnodesit = globalcollectednodes->begin(); gnodesit != globalcollectednodes->end(); ++gnodesit)
        {
          const int nodegid = (*gnodesit);
          const int nodelid = fluiddis->NodeRowMap()->LID(nodegid);
          if (nodelid >= 0)
            allcollectednodes->insert(nodegid);
        }
      }
    } // loop over layers of elements

    //-----------------------------------------------------------------------------------------------
    // If a node does not have 8 elements in the allcollected elements it must be a the surface
    // and therefore gets added to the surfacenodes set. This set is redundantly available and
    // mereley knows a node's position and velocities
    //-----------------------------------------------------------------------------------------------
    Teuchos::RCP<vector<LINALG::Matrix<3,2> > > surfacenodes = rcp(new vector<LINALG::Matrix<3,2> >);

    set<int>::const_iterator nodeit;
    for (nodeit = allcollectednodes->begin(); nodeit != allcollectednodes->end(); ++nodeit)
    {
      const int nodelid   = fluiddis->NodeRowMap()->LID(*nodeit);
      DRT::Node* node     = fluiddis->lRowNode(nodelid);
      DRT::Element** eles = node->Elements();
      int elementcount    = 0;
      for (int iele = 0; iele < node->NumElement(); ++iele)
      {
        DRT::Element* ele = eles[iele];
        set<int>::const_iterator foundit = allcollectedelements->find(ele->Id());
        if (foundit != allcollectedelements->end())
          elementcount++;
      }

      if (elementcount < 8)
      {
        LINALG::Matrix<3,2> coordandvel;
        const double* coord = node->X();
        const std::vector<int> dofids = (*dofset).Dof(node);
        for (int i = 0; i < 3; ++i)
        {
          const int lid = convel->Map().LID(dofids[i]);
          if (lid < 0)
            dserror("Dof %d not found on this proc.", dofids[i]);
          coordandvel(i,0) = coord[i];
          coordandvel(i,1) = (*convel)[lid];
        }
        surfacenodes->push_back(coordandvel);
      }
    }

    // Now the surfacenodes must be gathered to all procs
    {
      Teuchos::RCP<vector<LINALG::Matrix<3,2> > > mysurfacenodes = surfacenodes;
      surfacenodes = rcp(new vector<LINALG::Matrix<3,2> >);

      LINALG::Gather<LINALG::Matrix<3,2> >(*mysurfacenodes,*surfacenodes,numproc,allproc,fluiddis->Comm());
    }

    //----------------------------------------------------------------------------------------------
    // Here we manipulate the velocity vector. If a node is not in allnodes we find the nearest node
    // in surface nodes and use its velocity instead.
    //----------------------------------------------------------------------------------------------
    for(int lnodeid = 0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      DRT::Node* lnode = fluiddis->lRowNode(lnodeid);

      LINALG::Matrix<3,1> fluidvel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);

      // extract velocity values (no pressure!) from global velocity vector
      for (int i = 0; i < 3; ++i)
      {
        const int lid = convel->Map().LID(dofids[i]);
        if (lid < 0)
          dserror("Dof %d not found on this proc.");
        fluidvel(i) = (*convel)[lid];
      }

      set<int>::const_iterator foundit = allcollectednodes->find(lnode->Id());
      if (foundit == allcollectednodes->end())
      {
        // find closest node in surfacenodes
        LINALG::Matrix<3,2> closestnodedata(true);
        {
          LINALG::Matrix<3,1> nodecoord;
          const double* coord = lnode->X();
          for (int i = 0; i < 3; ++i)
            nodecoord(i) = coord[i];
          double mindist = 1.0e19;

          //--------------------------------
          // due to the PBCs the node might actually be closer to the
          // surfacenodes then would be calculated if one only considered
          // the actual position of the node. In order to find the
          // smallest distance the node is copied along all PBC directions
          //
          //   +------------------+ - - - - - - - - - -+
          //   +             II   +
          //   +   x        I  I  +    y               +
          //   +             II   +
          //   +------------------+ - - - - - - - - - -+
          //         original           copy
          //
          //   x: current node
          //   y: copy of current node
          //   I: interface
          //   +: pbc

          if (planenormal.size() > 3)
            dserror("Sorry, but currently a maximum of three periodic boundary conditions are supported by the combustion reinitializer.");

          // since there is no stl pow(INT, INT) function, we calculate it manually
          size_t looplimit = 1;
          for (size_t i = 0; i < planenormal.size(); ++i)
            looplimit *= 2;

          for (size_t ipbc = 0; ipbc < looplimit; ++ipbc)
          {
            LINALG::Matrix<3,1> tmpcoord(nodecoord);

            // determine which pbcs have to be applied
            //
            // loopcounter | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
            // ------------+---+---+---+---+---+---+---+---+
            //  first PBC  |     x       x       x       x
            // second PBC  |         x   x           x   x
            //  third PBC  |                 x   x   x   x
            //
            // this is equivalent to the binary representation
            // of the size_t
            if (ipbc & 0x01)
            {
              const double pbclength = globalmaxs[0] - globalmins[0];
              if (nodecoord(planenormal[0]) > globalmins[0] + pbclength/2.0)
                tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) - pbclength;
              else
                tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) + pbclength;
            }
            if (ipbc & 0x02)
            {
              const double pbclength = globalmaxs[1] - globalmins[1];
              if (nodecoord(planenormal[1]) > globalmins[1] + pbclength/2.0)
                tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) - pbclength;
              else
                tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) + pbclength;
            }
            if (ipbc & 0x04)
            {
              const double pbclength = globalmaxs[2] - globalmins[2];
              if (nodecoord(planenormal[2]) > globalmins[2] + pbclength/2.0)
                tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) - pbclength;
              else
                tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) + pbclength;
            }

            for (size_t i = 0; i < surfacenodes->size(); ++i)
            {
              const double dist = sqrt((tmpcoord(0)-(*surfacenodes)[i](0,0))*(tmpcoord(0)-(*surfacenodes)[i](0,0)) + (tmpcoord(1)-(*surfacenodes)[i](1,0))*(tmpcoord(1)-(*surfacenodes)[i](1,0)) + (tmpcoord(2)-(*surfacenodes)[i](2,0))*(tmpcoord(2)-(*surfacenodes)[i](2,0)));
              if (dist < mindist)
              {
                mindist = dist;
                closestnodedata = (*surfacenodes)[i];
              }
            }
          }
        }

        // write new velocities to the current node's dofs
        for (int icomp=0; icomp<3; ++icomp)
        {
          int lid = convel->Map().LID(dofids[icomp]);
          if (lid < 0)
            dserror("Dof %d not found on this proc.");
          const int err = conveltmp->ReplaceMyValues(1,&closestnodedata(icomp,1),&lid);
          if (err)
            dserror("could not replace values for convective velocity");
        }
      }
      else
      {
        for (int icomp=0; icomp<3; ++icomp)
        {
          int lid = convel->Map().LID(dofids[icomp]);
          if (lid < 0)
            dserror("Dof %d not found on this proc.");
          const int err = conveltmp->ReplaceMyValues(1,&fluidvel(icomp),&lid);
          if (err)
            dserror("could not replace values for convective velocity");
        }
      }
    }

    return conveltmp;
  }
  else
    return convel;
}


/* -------------------------------------------------------------------------------*
 * Restart (fluid is solved before g-func)                               rasthofer|
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Restart(int step)
{
  if (Comm().MyPID()==0)
    std::cout << "Restart of combustion problem" << std::endl;

  // restart of scalar transport (G-function) field
  ScaTraField().ReadRestart(step);

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  //--------------------------
  // write output to Gmsh file
  //--------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_after_restart", Step(), 701, true, gfuncdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phin \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phin(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinm \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinm(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Convective Velocity \" {" << endl;
    // draw vector field 'Convective Velocity' for every element
    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField().ConVel(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  gmshfilecontent.close();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  Teuchos::RCP<COMBUST::FlameFront> flamefrontOld = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis,ScaTraField().PBCmap()));

  // export phi n-1 vector from scatra dof row map to fluid node column map
  const Teuchos::RCP<Epetra_Vector> phinrow = rcp(new Epetra_Vector(*fluiddis->NodeRowMap()));
  if (phinrow->MyLength() != ScaTraField().Phin()->MyLength())
    dserror("vectors phinrow and phin must have the same length");
  *phinrow = *ScaTraField().Phin();
  const Teuchos::RCP<Epetra_Vector> phincol = rcp(new Epetra_Vector(*fluiddis->NodeColMap()));
  LINALG::Export(*phinrow,*phincol);

  // reconstruct old flame front
  //flamefrontOld->ProcessFlameFront(ScaTraField().Phin());
  flamefrontOld->ProcessFlameFront(phincol);

  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle =
      rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> dummy =
        rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  interfacehandle->UpdateInterfaceHandle();
  dummy->UpdateInterfaceHandle();
  FluidField().ImportInterface(interfacehandle,dummy);

  // restart of fluid field
  FluidField().ReadRestart(step);

  // reset interface for restart
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  interfacehandle_->UpdateInterfaceHandle();
  dummy->UpdateInterfaceHandle();
  //-------------------
  // write fluid output
  //-------------------
  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  FluidField().Output();
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  SetTimeStep(FluidField().Time(),step);

  UpdateTimeStep();

  return;
}

/* -------------------------------------------------------------------------------*
 * Restart (g-func is solved before fluid)                               rasthofer|
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::RestartNew(int step)
{
  if (Comm().MyPID()==0)
    std::cout << "Restart of combustion problem" << std::endl;

#ifdef BCF_IGNORE_SCATRA_RESTART
  Teuchos::RCP<Epetra_Vector> oldphinp = rcp(new Epetra_Vector(*(ScaTraField().Phinp())));
  Teuchos::RCP<Epetra_Vector> oldphin = rcp(new Epetra_Vector(*(ScaTraField().Phin())));
#endif

  // restart of scalar transport (G-function) field
  ScaTraField().ReadRestart(step);

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  //--------------------------
  // write output to Gmsh file
  //--------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_after_restart", Step(), 701, true, gfuncdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phin \" {" << endl;
    // draw scalar field 'Phin' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phin(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinm \" {" << endl;
    // draw scalar field 'Phinm' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinm(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "Convective Velocity \" {" << endl;
//    // draw vector field 'Convective Velocity' for every element
//    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField().ConVel(),gmshfilecontent);
//    gmshfilecontent << "};" << endl;
  }
  gmshfilecontent.close();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  //Teuchos::RCP<COMBUST::FlameFront> flamefrontOld = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis,));

//  // export phi n vector from scatra dof row map to fluid node column map
//  const Teuchos::RCP<Epetra_Vector> phinrow = rcp(new Epetra_Vector(*fluiddis->NodeRowMap()));
//  if (phinrow->MyLength() != ScaTraField().Phin()->MyLength())
//    dserror("vectors phinrow and phin must have the same length");
//  *phinrow = *ScaTraField().Phin();
//  const Teuchos::RCP<Epetra_Vector> phincol = rcp(new Epetra_Vector(*fluiddis->NodeColMap()));
//  LINALG::Export(*phinrow,*phincol);

  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  interfacehandle_->UpdateInterfaceHandle();

  // reconstruct old flame front
  //flamefrontOld->ProcessFlameFront(phincol);

  // build interfacehandle using old flame front
  // TODO @Martin Test + Kommentar
  // remark: interfacehandleN = interfacehandleNP, weil noch aeltere Information nicht vorhanden

  // TODO remove old code when new code tested
  //Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOld =
  //  rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  //interfacehandleOld->UpdateInterfaceHandle();
  //FluidField().ImportInterface(interfacehandleOld);

//  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle =
//      rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));//flamefrontOld));
//  Teuchos::RCP<COMBUST::InterfaceHandleCombust> dummy =
//        rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));//flamefrontOld));
//  interfacehandle->UpdateInterfaceHandle();
//  dummy->UpdateInterfaceHandle(); // dummy because second function argument required

    // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  //FluidField().ImportInterface(interfacehandle,dummy);
  FluidField().ImportInterface(interfacehandle_,interfacehandle_);
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  // restart of fluid field
  FluidField().ReadRestart(step);

#ifdef BCF_IGNORE_SCATRA_RESTART
  // now overwrite restart phis w/ the old phis
  ScaTraField().Phinp()->Update(1.0, *(oldphinp), 0.0);
  ScaTraField().Phin()->Update(1.0, *(oldphin), 0.0);
  ScaTraField().ComputeTimeDerivative();
  ScaTraField().Phidtn()->Update(1.0,*(ScaTraField().Phidtnp()),0.0);
  ScaTraField().Phinm()->Update(1.0,*(ScaTraField().Phin()),0.0);

  // additionally we need to update the interfacehandle and flamefront
  // or later on the computeVolume function will return the old volume
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());
  interfacehandle_->UpdateInterfaceHandle();
  FluidField().ImportFlameFront(flamefront_);
  FluidField().ImportInterface(interfacehandle_,interfacehandle_);
  FluidField().ImportFlameFront(Teuchos::null);
#endif

  //-------------------
  // write fluid output
  //-------------------
//  interfacehandle_->UpdateInterfaceHandle();
//  interfacehandleold_->UpdateInterfaceHandle();
//  // show flame front to fluid time integration scheme
//  FluidField().ImportFlameFront(flamefront_);
//  FluidField().Output();
//  // delete fluid's memory of flame front; it should never have seen it in the first place!
//  FluidField().ImportFlameFront(Teuchos::null);

  SetTimeStep(FluidField().Time(),step);

//  UpdateTimeStep();

  return;
}

#endif // #ifdef CCADISCRET
