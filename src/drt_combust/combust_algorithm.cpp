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

#include "combust_algorithm.H"
#include "combust_defines.H"
#include "combust_flamefront.H"
#include "two_phase_defines.H"
#include "combust_reinitializer.H"
#include "combust_fluidimplicitintegration.H"
#include "combust_reinitialization_pde.H"
#include "../drt_fluid/turbulence_statistic_manager.H"
#include "../drt_fluid/turbulence_statistics_mean_general.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_io/io_gmsh.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
//#include "../drt_lib/standardtypes_cpp.H" // required to use mathematical constants such as PI


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
COMBUST::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn, const Teuchos::ParameterList& solverparams)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn,false,"scatra",solverparams),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  convtol_(combustdyn.get<double>("CONVTOL")),
  stepbeforereinit_(false),
  stepreinit_(false),
  combusttype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  reinitaction_(DRT::INPUT::IntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
//  reinitaction_(combustdyn.sublist("COMBUSTION GFUNCTION").get<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION")),
  reinitinitial_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINIT_INITIAL")),
  reinitinterval_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("REINITINTERVAL")),
  reinitband_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITBAND")),
  reinitbandwidth_(combustdyn.sublist("COMBUSTION GFUNCTION").get<double>("REINITBANDWIDTH")),
  reinit_output_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINIT_OUTPUT")),
  volcorrection_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITVOLCORRECTION")),
  evaltimeratio_(1.0),
  extract_interface_vel_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"EXTRACT_INTERFACE_VEL")),
  convel_layers_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("NUM_CONVEL_LAYERS")),
  gmshoutput_(combustdyn.get<bool>("GMSH_OUTPUT")),
  combustdyn_(combustdyn),
  flamefront_(Teuchos::null),
  reinit_pde_(Teuchos::null),
  transport_vel_(DRT::INPUT::IntegralValue<INPAR::COMBUST::TransportVel>(combustdyn.sublist("COMBUSTION GFUNCTION"),"TRANSPORT_VEL")),
  transport_vel_no_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("TRANSPORT_VEL_FUNC")),
  restart_(false)
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

  // TODO: scatra and fluid/combustion time-integration methods do not have to fit, or shall they? - winklmaier
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


  // FGI vectors are initialized
  velnpi_ = rcp(new Epetra_Vector(FluidField().StdVelnp()->Map()),true);//*fluiddis->DofRowMap()),true);
  velnpi_->Update(1.0,*FluidField().StdVelnp(),0.0);
  phinpi_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
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
  flamefront_->UpdateOldInterfaceHandle();

  volume_start_ = ComputeVolume();

  if(reinitinitial_)
  {
    //--------------------------------------------------------------
    // initial reinitialization by geometric distance computation
    // remark: guarantee (periodic) initial signed distance property
    //--------------------------------------------------------------
    if (reinitaction_ != INPAR::COMBUST::reinitaction_none)
    {
      //cou t<< "reinitialization at timestep 0 switched off!" << endl;

    // get my flame front (boundary integration cells)
    std::map<int, GEO::BoundaryIntCells> myflamefront = flamefront_->InterfaceHandle()->BoundaryIntCells();
#ifdef PARALLEL
      // export flame front (boundary integration cells) to all processors
      flamefront_->ExportFlameFront(myflamefront);
#endif
      // reinitialize initial level set field by geometric distance computation
      COMBUST::Reinitializer reinitializer(
          combustdyn_,
          ScaTraField(),
          myflamefront,
          ScaTraField().Phinp(),
          true);

      if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary)
      {
        // reset phin vector in ScaTra time integration scheme to phinp vector
        *ScaTraField().Phin() = *ScaTraField().Phinp();
      }

    // update flame front according to reinitialized G-function field
    flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

      // pointer not needed any more
      stepreinit_ = false;
    }
  }
  if(reinitaction_  == INPAR::COMBUST::reinitaction_sussman )
  {
    reinit_pde_ = rcp(new COMBUST::ReinitializationPDE());
  }

  //---------------------------------------------------
  // set initial fluid field (based on standard dofset)
  //---------------------------------------------------
  // update fluid interface with flamefront
  FluidField().ImportFlameFront(flamefront_,false);
  // read parameters for initial field
  const INPAR::COMBUST::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(
      combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD");
  const int initfuncno = combustdyn_.sublist("COMBUSTION FLUID").get<int>("INITFUNCNO");
  // set initial flow field based on standard dofs only
  FluidField().SetInitialFlowField(initfield, initfuncno);
  // clear fluid's memory to flamefront
  FluidField().ImportFlameFront(Teuchos::null,false);

  //--------------------------------------------------------
  // incorporate the XFEM dofs into the fluid
  // remark: this includes setting the initial enriched dofs
  //--------------------------------------------------------
  // update fluid interface with flamefront
  FluidField().ImportFlameFront(flamefront_,true);
  // output fluid initial state
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary)
    FluidField().Output();
  // clear fluid's memory to flamefront
  FluidField().ImportFlameFront(Teuchos::null,true);

  if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion)
  {
    // extract convection velocity from fluid solution
    const Teuchos::RCP<const Epetra_Vector> convel = (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
                                               ?(FluidField().StdVelaf())
                                               :(FluidField().StdVelnp());

    if(transport_vel_ == INPAR::COMBUST::transport_vel_flamevel)
    {
      ScaTraField().SetVelocityField(ComputeFlameVel(convel,FluidField().DofSet()),
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    else if(transport_vel_ == INPAR::COMBUST::transport_vel_byfunct)
    {
      ScaTraField().SetVelocityField(OverwriteFluidVel(),
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    else dserror("no valid transport velocity for combust problems!");
  }

  // output G-function initial state
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary and
      DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false )
    ScaTraField().Output();

//  if (Comm().MyPID()==0)
//  {
//    std::string outfile = "/home/henke/simulations/results/study_timeint/2d_64_neumann_visc-3_theta05/radius";
//    //outfile.append(".oracles.-5h.flow_statistics");
//    std::ofstream title(outfile.c_str(),ios::out);
//    title << "# radius of circular flame\n";
//    title << "# reinitialized radius  original radius\n";
//    title << "# initial radius 0.025";
//    title << "\n";
//    title.flush();
//  }
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
  // however, calculate it if and only if the problem has not been restarted
  if(DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == true and restart_==false)
    SolveInitialStationaryProblem();

  // time loop
  while (NotFinished())
  {
    // redistribute
    Redistribute();

    // prepare next time step; update field vectors
    PrepareTimeStep();

    // Fluid-G-function-Interaction loop
    // remark: - the G-function is solved first, then the fluid is solved
    //         - this guarantees that the interface geometry and the fluid solution match
    //         - this is a prerequisite for the semi-Lagrangean time integration method
    //         - changing this order would affect the restart
    //         -> do not change this order!
    while (NotConvergedFGI())
    {
      // prepare Fluid-G-function iteration
      PrepareFGIteration();

      // solve linear G-function equation
      //cout << "---------------------" << endl;
      //cout << "Level set deactivated!" << endl;
      //cout << "---------------------" << endl;
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
  volume_start_ = ComputeVolume();

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
  // unnecessary for stationary problems
  //UpdateInterface();

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
  printMassConservationCheck(volume_start_, volume_end);

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                schott 05/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoReinitialization()
{

  // compute current volume of minus domain
  const double volume_current_before = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start_, volume_current_before);


  // reinitialize Gfunc
  switch(reinitaction_)
  {
    case INPAR::COMBUST::reinitaction_sussman:
    {
      reinitialization_accepted_ = reinit_pde_->CallReinitialization(ScaTraField().Phinp());
      break;
    }
    case INPAR::COMBUST::reinitaction_signeddistancefunction:
    case INPAR::COMBUST::reinitaction_fastsigneddistancefunction:
      if(Comm().MyPID()==0) std::cout << "COMBUST::Algorithm: reinitialization via calculating signed distance functions" << std::endl;
      ReinitializeGfuncSignedDistance();
      reinitialization_accepted_ = true;
      break;
    case INPAR::COMBUST::reinitaction_none:
      if(Comm().MyPID()==0) std::cout << "No reinitialization chosen" << std::flush;
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
    // reinitilization is done every fgi -> this is neither good nor cheap,
    // but we didn't find a better solution...
    if (fgitermax_>1)
      cout << "WARNING: reinitialization done every FLI" << endl;

    // update flame front according to reinitialized G-function field
    flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phinp());
  }

  // TODO: the step used for output is switched for one
  if(reinit_output_ and reinitialization_accepted_ and !(reinitaction_ == INPAR::COMBUST::reinitaction_signeddistancefunction or reinitaction_ == INPAR::COMBUST::reinitaction_fastsigneddistancefunction))
  {
    reinit_pde_->OutputReinit(ScaTraField().Step(), ScaTraField().Time());
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
    std::map<int, GEO::BoundaryIntCells> myflamefront = flamefront_->InterfaceHandle()->BoundaryIntCells();

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
      flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phiaf());
    else
      flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phinp());


    // get my flame front (boundary integration cells)
    // remark: we generate a copy of the flame front here, which is not neccessary in the serial case
    std::map<int, GEO::BoundaryIntCells> myflamefront = flamefront_->InterfaceHandle()->BoundaryIntCells();

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
      flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phiaf());
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
      flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phinp());
    }
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

  if(transport_vel_no_ == -1) dserror("please set a function number for interface transport velocity!");

  // get fluid (Navier-Stokes) velocity vector in standard FEM configuration (no XFEM dofs)
  const Teuchos::RCP<Epetra_Vector> convel = FluidField().StdVelnp();
  // velocity function number = 1 (FUNCT1)
  const int velfuncno = transport_vel_no_;

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
      double value = DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(icomp,lnode->X(),FluidField().Time(),NULL);

      // scaling for stretching fluid example // schott
      //cout <<"Scaling with time-curve!!!!" << endl;
      //value *= cos(PI*FluidField().Time()/500.0);
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
  const Teuchos::RCP<Epetra_Vector> & convel = FluidField().StdVelnp();

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
      if (gmshoutput_)
      {
        LINALG::SerialDenseMatrix xyz(3,1);
        xyz(0,0) = lnode->X()[0];
        xyz(1,0) = lnode->X()[1];
        xyz(2,0) = lnode->X()[2];

        IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, nvec, xyz, gmshfilecontent);
      }
#endif

      //------------------------
      //      // get material parameters
      //      //------------------------
      //      // get list of adjacent elements of this node
      //      DRT::Element** elelist = lnode->Elements();
      //
      //      // get material from first (arbitrary!) element adjacent to this node
      //      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
      //      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
      //      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());
      //
      //      // density burnt domain
      //      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
      //      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      //      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
      //      const double rhoplus = matplus->Density();
      //
      //      // density unburnt domain
      //      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
      //      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      //      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
      //      const double rhominus = matminus->Density();
      //
      //      // laminar flame speed
      //      const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
      //      //---------------------------------------------
      //      // compute relative flame velocity at this node
      //      //---------------------------------------------
      //      // get phi value for this node
      //      const double gfuncval = (*phinp)[nodelid];
      //
      //      double speedfac = 0.0;
      //      if (gfuncval >= 0.0){ // interface or burnt domain -> burnt material
      //        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
      //        speedfac = sl * rhominus/rhoplus;
      //      }
      //      else{ // unburnt domain -> unburnt material
      //        // flame speed factor = laminar flame speed
      //        speedfac = sl;
      //      }


      //      cout << "WARNING: the velocity vector is set to zero at the moment" << endl;
      //      nvec.Clear();

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
    const Teuchos::RCP<const Epetra_Vector>& convel,
    const Teuchos::RCP<const DRT::DofSet>& dofset
    //const Teuchos::RCP<const Epetra_Map >& dbcmap
    )
{
  if (Comm().MyPID()==0)
  {
    if((DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false) and
       (DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD") == INPAR::COMBUST::initfield_zero_field))
      cout << "/!\\ Compute an initial stationary fluid solution to avoid a non-zero initial flame velocity" << endl;
  }

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
  const double marksteinlength  = combustdyn_.sublist("COMBUSTION FLUID").get<double>("MARKSTEIN_LENGTH");

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

#ifdef ORACLES
      // skip the part of the domain left of the expansion
      if (lnode->X()[0] < 0.0)
        continue;
#endif
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
      if (fabs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
      {
        // length of smoothed normal is zero at this node -> node must be on a singularity of the
        // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
        // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
        //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
        std::cout << "\n/!\\ phi gradient too small at node " << gid << " -> flame velocity is only the convective velocity" << std::endl;
        nvec.PutScalar(0.0);
      }
      else
      {
        nvec.Scale(-1.0/ngradphi);
      }

#ifdef COLLAPSE_FLAME
      nvec(0) = lnode->X()[0];
      nvec(1) = lnode->X()[1];
      nvec(2) = lnode->X()[2];
#ifdef COMBUST_2D
      nvec(2) = 0.0;
#endif
      const double norm = nvec.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0)
      {
        nvec.PutScalar(0.0);
        //dserror("norm of normal vector is zero!");
      }
      else
      {
        nvec.Scale(1.0/norm);
      }
#endif

#ifdef COMBUST_GMSH_NORMALFIELD
      if (gmshoutput_)
      {
        LINALG::SerialDenseMatrix xyz(3,1);
        xyz(0,0) = lnode->X()[0];
        xyz(1,0) = lnode->X()[1];
        xyz(2,0) = lnode->X()[2];

        IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, nvec, xyz, gmshfilecontent);
      }
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
      const double curv = (*curvature)[nodelid];

      double speedfac = 0.0;

#ifdef COMBUST_BUNSENBURNER
      //      if(((lnode->X()[1]>0.5 && lnode->X()[1]<=0.55) && (lnode->X()[0]>=0.45)) || ((lnode->X()[1]>0.5 && lnode->X()[1]<=0.55) && (lnode->X()[0]<=-0.45)))
      //      {
      //          if(lnode->X()[0]>0)
      //              speedfac = 10.0*sqrt((lnode->X()[1]-0.5)*(lnode->X()[1]-0.5)+(lnode->X()[0]-0.5)*(lnode->X()[0]-0.5));
      //          else
      //              speedfac = 10.0*sqrt((lnode->X()[1]-0.5)*(lnode->X()[1]-0.5)+(lnode->X()[0]+0.5)*(lnode->X()[0]+0.5));
      //      }
      if(lnode->X()[1]<=0.5 && lnode->X()[0]<0.51 && lnode->X()[0]>-0.51)
      {
        speedfac = 0.0;
      } //2D
      //      if(lnode->X()[1]<=0.0 && lnode->X()[0]<0.51 && lnode->X()[0]>-0.51)
      //      {
      //          speedfac = 0.0;
      //      } //2D einfach
      //      if(lnode->X()[2]<=0.5 && sqrt(lnode->X()[0]*lnode->X()[0]+lnode->X()[1]*lnode->X()[1])<0.51)
      //      {
      //          speedfac = 0.0;
      //      }  //3D
#endif //BUNSENBURNER

      double wallfac = 1.0;
#ifdef ORACLES
      //---------------------------------------------------------
      // blend the flame speed close to walls for ORACLES problem
      //---------------------------------------------------------
      //    wall
      // 1.0 |     ______
      //     |    /
      // 0.0 |___/ ,
      //     |     H/6
      const double wallzone = 0.0299/6.0;
      if (lnode->X()[0] > 0.0) // inside combustion chamber
      {
        if ( (0.0653-abs(lnode->X()[1])) < wallzone or // close to top or bottom wall
                         lnode->X()[0]   < wallzone )  // close to step
        {
          // wall factor is 0 at the wall and 1 at H/6 or further away from the wall
          wallfac = 6.0/0.0299 * std::min(0.0653-abs(lnode->X()[1]),lnode->X()[0]);
        }
      }
#endif

      if (XFEM::plusDomain(gfuncval) == true){ // interface or burnt domain -> burnt material
        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
        speedfac = sl*(1.0-marksteinlength*curv)* rhominus/rhoplus;
      }
      else{ // unburnt domain -> unburnt material
        // flame speed factor = laminar flame speed
        speedfac = sl*(1.0-marksteinlength*curv);
      }

      LINALG::Matrix<3,1> flvelrel(true);
      for (int icomp=0; icomp<3; ++icomp)
        flvelrel(icomp) = wallfac * speedfac * nvec(icomp);

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

#ifdef ORACLES
      // set relative convective velocity to zero at the walls
      // remark: - the wall factor does the same thing, this is just to be sure
      //         - this is a hack for ORACLES
      //         - this guarantees a zero transport velocity for no-slip Dirichlet walls
      //         - this could be done in clean way by using a vector holding the fluid Dirichlet conditions
      if ( (abs(lnode->X()[1])-0.0653) < 1.0E-9 or // top and bottom wall combustion chamber
          ( abs(lnode->X()[0]) < 1.0E-9 and abs(lnode->X()[1]) > (0.005-1.0E-9) ) ) // walls above and below step
        flvelrel.Clear();

      // this is an attempt to do this in a clean general way; it failed
      // set Dirichlet BC
      //for (int icomp=0; icomp<3; ++icomp)
      //{
      //  const int gid = dofids[icomp];
      //  cout << gid << endl;
      //  if(dbcmap->MyGID(gid))
      //  {
      //    cout << "DBC overwritten" << endl;
      //    cout << nodecoord(0,0) << endl;
      //    cout << nodecoord(1,0) << endl;
      //    cout << nodecoord(2,0) << endl;
      //    flvelrel(icomp) = 0.0;
      //  }
      //}
#endif

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
  bool notconverged = true;

  if (fgitermax_ == 0)
    dserror("A maximum of 0 FGI is not sensible!!!");

  if (fgiter_ > 0) // nothing to do for FGI = 0
  {
    const Teuchos::RCP<const Epetra_Vector> velnpip = FluidField().StdVelnp();
    const Teuchos::RCP<const Epetra_Vector> phinpip = ScaTraField().Phinp();

    double velnormL2 = 1.0;
    double gfuncnormL2 = 1.0;

    velnpip->Norm2(&velnormL2);
    phinpip->Norm2(&gfuncnormL2);

    if (velnormL2 < 1e-5) velnormL2 = 1.0;
    if (gfuncnormL2 < 1e-5) gfuncnormL2 = 1.0;

    double fgvelnormL2 = 1.0;
    double fggfuncnormL2 = 1.0;

    // compute increment and L2-norm of increment
    RCP<Epetra_Vector> incvel = rcp(new Epetra_Vector(velnpip->Map()),true);
    incvel->Update(1.0,*velnpip,-1.0,*velnpi_,0.0);
    incvel->Norm2(&fgvelnormL2);

    RCP<Epetra_Vector> incgfunc = rcp(new Epetra_Vector(phinpip->Map(),true));//*ScaTraField().Discretization()->DofRowMap()),true);
    incgfunc->Update(1.0,*phinpip,-1.0,*phinpi_,0.0);
    incgfunc->Norm2(&fggfuncnormL2);

    if (fgitermax_ > 1)
    {
      if (Comm().MyPID()==0)
      {
        printf("\n|+------------------------ FGI ------------------------+|");
        printf("\n|iter/itermax|----tol-[Norm]--|-fluid inc--|-g-func inc-|");
        printf("\n|   %2d/%2d    | %10.3E[L2] | %10.3E | %10.3E |",
            fgiter_,fgitermax_,convtol_,fgvelnormL2/velnormL2,fggfuncnormL2/gfuncnormL2);
        printf("\n|+-----------------------------------------------------+|\n");

        if (fgiter_ == fgitermax_)
        {
          printf("|+---------------- not converged ----------------------+|");
          printf("\n|+-----------------------------------------------------+|\n");
        }
      } // end if processor 0 for output

      if ((fgvelnormL2/velnormL2 <= convtol_) and (fggfuncnormL2/gfuncnormL2 <= convtol_))
        notconverged = false;
    } // end if fgitermax > 1

    if (fgiter_ >= fgitermax_)
      notconverged = false;

    if (!notconverged)
      FluidField().ClearTimeInt(); // clear XFEM time integration

    // update fgi-vectors
    velnpi_->Update(1.0,*velnpip,0.0);
    phinpi_->Update(1.0,*phinpip,0.0);
  }

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
      FluidField().StdVelnp(),
      Teuchos::null,
      Teuchos::null,
      Teuchos::null,
      FluidField().DofSet(),
      FluidField().Discretization()
    );

    // Transfer history vector only for subgrid-velocity
    //ScaTraField().SetVelocityField(
    //    FluidField().ExtractInterfaceVeln(),
    //    FluidField().Hist(),
    //    Teuchos::null,
    //    FluidField().DofSet(),
    //    FluidField().Discretization()
    //);
    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added

    // extract convection velocity from fluid solution
    const Teuchos::RCP<const Epetra_Vector> convel = FluidField().StdVelnp();

    ScaTraField().SetVelocityField(
//        OverwriteFluidVel(),
        //FluidField().ExtractInterfaceVeln(),
        ComputeFlameVel(convel,FluidField().DofSet()),
        Teuchos::null,
        Teuchos::null,
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
  if (Comm().MyPID()==0 and fgitermax_>1)
  {
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

  // update fluid interface with flamefront
  FluidField().ImportFlameFront(flamefront_,true);

  // solve nonlinear Navier-Stokes equations
  FluidField().NonlinearSolve();

  // clear fluid's memory to flamefront
  FluidField().ImportFlameFront(Teuchos::null,false);


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

  // get the convel at the correct time
  const Teuchos::RCP<const Epetra_Vector> convel = (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
                                                   ?(FluidField().StdVelaf())
                                                   :(FluidField().StdVelnp());
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
      ScaTraField().SetVelocityField(ManipulateFluidFieldForGfunc(convel, FluidField().DofSet()),
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    else // transfer the computed velocity as convective velocity to the g-function solver
    {
      ScaTraField().SetVelocityField(convel,
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());

      // Transfer history vector only for subgrid-velocity
      //ScaTraField().SetVelocityField(
      //    FluidField().ExtractInterfaceVeln(),
      //    FluidField().Hist(),
      //    Teuchos::null,
      //    Teuchos::null,
      //    FluidField().DofSet(),
      //    FluidField().Discretization()
      //);
    }
    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added

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

    if(transport_vel_ == INPAR::COMBUST::transport_vel_flamevel)
    {
      ScaTraField().SetVelocityField(ComputeFlameVel(convel,FluidField().DofSet()),
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    else if(transport_vel_ == INPAR::COMBUST::transport_vel_byfunct)
    {
      ScaTraField().SetVelocityField(OverwriteFluidVel(),
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    else dserror("no valid transport velocity for combust problems!");
    break;
  }
  default:
    dserror("unknown type of combustion problem");
    break;
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
  //overwrite old interfacehandle before updating flamefront in first FGI
  if (fgiter_<=1)
    flamefront_->UpdateOldInterfaceHandle();

  // update flame front according to evolved G-function field
  // remark: for only one FGI iteration, 'phinpip_' == ScaTraField().Phin()
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") == INPAR::FLUID::timeint_afgenalpha)
    flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phiaf());
  else
    flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField().Phinp());

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
    if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      ScaTraField().SetVelocityField(FluidField().StdVeln(), //as the fluid has already been updated, we need veln here
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     FluidField().DofSet(),
                                     FluidField().Discretization());
    }
    ScaTraField().UpdateReinit();
    // for genalpha there is no need to write velaf into the scatra field, as it will be done at the beginning of the next timestep anyways
  }
  else
  {
    //if(Comm().MyPID()==0)
    //  cout << "Update" << endl;
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
  FluidField().ImportFlameFront(flamefront_,false);
  // write fluid output
  FluidField().Output();
  // clear fluid's memory to flamefront
  FluidField().ImportFlameFront(Teuchos::null,false);


  // causes error in DEBUG mode (trueresidual_ is null)
  //FluidField().LiftDrag();
  ScaTraField().Output();

  // we have to call the output of averaged fields for scatra separately
  if (FluidField().TurbulenceStatisticManager() != Teuchos::null)
    FluidField().TurbulenceStatisticManager()->DoOutputForScaTra(ScaTraField().DiscWriter(),ScaTraField().Step());

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
      std::cout << "           mass conservation"            << endl;
      std::cout << " initial mass: " << volume_start << endl;
      std::cout << " final mass:   " << volume_end   << endl;
      std::cout << " mass loss:    " << massloss << "%" << endl;
      std::cout << "---------------------------------------" << endl;
    }
    else
    {
      cout << " there is no 'minus domain'! -> division by zero checking mass conservation" << endl;
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
  double myvolume = flamefront_->InterfaceHandle()->ComputeVolumeMinus();

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
    std::cout << "Correcting volume of minus(-) domain by " << voldelta << " ... " << std::flush;

  // compute negative volume of discretization on this processor
  double myarea = flamefront_->InterfaceHandle()->ComputeSurface();
  double sumarea = 0.0;

  // sum volumes on all processors
  // remark: ifndef PARALLEL sumvolume = myvolume
  Comm().SumAll(&myarea,&sumarea,1);

  // This is a guess on how thick a layer needs to be added to the surface of the minus domain.
  // Due to $ \grad \phi \approx 1 $ this also happens to be the value that needs to be subtracted
  // of all phis. To make sure that \grad \phi really is close to 1 this function should only be
  // called after a reinitialization.
  const double thickness = -voldelta / sumarea;

  Teuchos::RCP<Epetra_Vector> phin = ScaTraField().Phin();

  Teuchos::RCP<Epetra_Vector> one = Teuchos::rcp(new Epetra_Vector(phin->Map()));
  one->PutScalar(1.0);

  if (phin != Teuchos::null)
    phin->Update(thickness, *one, 1.0);

  Teuchos::RCP<Epetra_Vector> phinp = ScaTraField().Phinp();
  if (phinp != Teuchos::null)
    phinp->Update(thickness, *one, 1.0);

  Teuchos::RCP<Epetra_Vector> phinm = ScaTraField().Phinm();
  if (phinm != Teuchos::null)
    phinm->Update(thickness, *one, 1.0);

  // REMARK:
  // after the reinitialization we update the flamefront in the usual sense
  // this means that we modify phi-values if necessary -> default boolian modifyPhiVectors = true

  // update flame front according to reinitialized G-function field
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  if (Comm().MyPID() == 0)
    std::cout << "done" << std::endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 |                                                                                rasthofer 08/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> COMBUST::Algorithm::ManipulateFluidFieldForGfunc(
    const Teuchos::RCP<const Epetra_Vector>& convel,
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

            if (XFEM::plusDomain((*phinp)[nodelid]) == false)
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
 | Restart a combustion problem                                          rasthofer|
 | remark: G-function is solved before fluid                                      |
 |         switching the order would affect the restart                           |
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Restart(int step, const bool restartscatrainput, const bool restartfromfluid)
{
  if (Comm().MyPID()==0)
  {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "| restart of combustion problem             |" << std::endl;
    if (restartscatrainput)
      std::cout << "| restart with scalar field from input file |" << endl;
    if (restartfromfluid)
      std::cout << "| restart from standard fluid problem       |" << endl;
    std::cout << "---------------------------------------------" << std::endl;
  }
  if (restartfromfluid and !restartscatrainput)
    dserror("scalar field must be read from input file for restart from standard fluid");

  // read level-set field from input file instead of restart file
  Teuchos::RCP<Epetra_Vector> oldphidtn  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphin    = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphinm   = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphinp   = Teuchos::null;
  if (restartscatrainput)
  {
    oldphin  = rcp(new Epetra_Vector(*(ScaTraField().Phin())));
    oldphinp = rcp(new Epetra_Vector(*(ScaTraField().Phinp())));
    if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
      oldphidtn= rcp(new Epetra_Vector(*(ScaTraField().Phidtn())));
    else
      oldphinm = rcp(new Epetra_Vector(*(ScaTraField().Phinm())));
  }

  // restart of scalar transport (G-function) field
  if (!restartfromfluid) // not if restart is done from standard fluid field; there is no scalar field
    ScaTraField().ReadRestart(step);

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());

  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_,true);
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null,false);

  // restart of fluid field
  FluidField().ReadRestart(step);

  // read level-set field from input file instead of restart file
  if (restartscatrainput)
  {
    if (Comm().MyPID()==0)
      std::cout << "---  overwriting scalar field with field from input file... " << std::flush;
    // now overwrite restart phis w/ the old phis
    ScaTraField().Phinp()->Update(1.0, *(oldphinp), 0.0);
    ScaTraField().Phin() ->Update(1.0, *(oldphin),  0.0);
    if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      ScaTraField().Phidtn() ->Update(1.0, *(oldphidtn),  0.0);
      ScaTraField().ComputeIntermediateValues();
    }
    else
    {
      ScaTraField().Phidtn()->PutScalar(0.0);
      ScaTraField().ComputeTimeDerivative();
      ScaTraField().Phidtn()->Update(1.0,*(ScaTraField().Phidtnp()),0.0);
      ScaTraField().Phinm() ->Update(1.0,*(ScaTraField().Phin()),   0.0);
    }
    if (Comm().MyPID()==0)
      std::cout << "done" << std::endl;

    //-------------------------------------------------------------
    // fill flamefront conforming to restart state
    //-------------------------------------------------------------
    flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());
    flamefront_->UpdateOldInterfaceHandle();

    // show flame front to fluid time integration scheme
    FluidField().ImportFlameFront(flamefront_,true);
    // delete fluid's memory of flame front; it should never have seen it in the first place!
    FluidField().ImportFlameFront(Teuchos::null,false);

    if (gmshoutput_)
    {
      //--------------------------
      // write output to Gmsh file
      //--------------------------
      const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_restart_before_reinit", Step(), 701, true, gfuncdis->Comm().MyPID());
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
      if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "Phidtn \" {" << endl;
        // draw scalar field 'Phidtn' for every element
        IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phidtn(),gmshfilecontent);
        gmshfilecontent << "};" << endl;
      }
      else
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
      std::cout << " done" << endl;
    }

    if (!restartfromfluid) // default: restart from an XFEM problem, not from a standard fluid problem
    {
      // get my flame front (boundary integration cells)
      std::map<int, GEO::BoundaryIntCells> myflamefront = flamefront_->InterfaceHandle()->BoundaryIntCells();
#ifdef PARALLEL
      // export flame front (boundary integration cells) to all processors
      flamefront_->ExportFlameFront(myflamefront);
#endif
      // reinitialize G-function (level set) field
      COMBUST::Reinitializer reinitializer(
          combustdyn_,
          ScaTraField(),
          myflamefront,
          ScaTraField().Phinp(),
          true);

      *ScaTraField().Phin() = *ScaTraField().Phinp();
      if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_one_step_theta)
        *ScaTraField().Phinm() = *ScaTraField().Phinp();
      //TODO: what should be done for gen alpha and does not phidtx have to be updated too?

      // update flame front according to reinitialized G-function field
      flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());
    }
  }

  phinpi_->Update(1.0,*ScaTraField().Phin(),0.0);

  if (gmshoutput_)
  {
    //--------------------------
    // write output to Gmsh file
    //--------------------------
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_restart_after_reinit", Step(), 701, true, gfuncdis->Comm().MyPID());
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
    if (ScaTraField().MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phidtn \" {" << endl;
      // draw scalar field 'Phidtn' for every element
      IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phidtn(),gmshfilecontent);
      gmshfilecontent << "};" << endl;
    }
    else
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
    std::cout << " done" << endl;
  }

  //FluidField().Output();
  //ScaTraField().Output();

  // set time in scalar transport time integration scheme
  if(restartfromfluid)
    ScaTraField().SetTimeStep(FluidField().Time(),step);

  SetTimeStep(FluidField().Time(),step);

  //UpdateTimeStep();

  restart_ = true;

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: redistribute scatra and fluid discretization                        rasthofer 07/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Redistribute()
{
  const double ratiolimitfac = combustdyn_.get<double>("PARALLEL_REDIST_RATIO_FAC");

  if (Comm().NumProc() > 1 and ratiolimitfac > 1.0)
  {
    // for ease of changing some constants are set here
    double normalnodeweight   = 10.0;
    double enrichednodeweight = 2000.0;
    double slavenodeweight    = 1.0;

    // determine whether a redistribution is necessary
    // by comparing min and max evaltime of all procs
    double myprocevaltime  = FluidField().EvalTime();
    double minprocevaltime = 0.0;
    double maxprocevaltime = 0.0;

    Comm().MinAll(&myprocevaltime, &minprocevaltime, 1);
    Comm().MaxAll(&myprocevaltime, &maxprocevaltime, 1);

    // Calculate the evaluation time ratio by which a further redistribution is determined.
    // For this the minimum ratio since the last redistribution is used
    if (maxprocevaltime / minprocevaltime < evaltimeratio_)
      evaltimeratio_ = maxprocevaltime / minprocevaltime;

    if (minprocevaltime <= 0.0)
    {
      if (Comm().MyPID() == 0)
        cout << "The max / min ratio could not be determined --> Not redistributing" << endl;
    }
    else if (maxprocevaltime / minprocevaltime < ratiolimitfac * evaltimeratio_)
    {
      if (Comm().MyPID() == 0)
        cout << "The max / min ratio is " << maxprocevaltime / minprocevaltime << " < " << ratiolimitfac * evaltimeratio_ << " --> Not redistributing" << endl;
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        cout << "-------------------------------Redistributing-------------------------------" << endl;
        cout << "The max / min ratio is " << maxprocevaltime / minprocevaltime << " > " << ratiolimitfac * evaltimeratio_ << " --> Redistributing" << endl;
      }
      //--------------------------------------------------------------------------------------
      // Building graph for later use by parmetis
      //--------------------------------------------------------------------------------------
      const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
      const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();
      const Epetra_Map* noderowmap = fluiddis->NodeRowMap();
      Teuchos::RCP<map<int,vector<int> > > allcoupledcolnodes = ScaTraField().PBC()->ReturnAllCoupledColNodes();

      // weights for graph partition
      Epetra_Vector weights(*noderowmap,false);
      weights.PutScalar(normalnodeweight);

      // loop elements and raise the weight of nodes of cut elements to 1000
      {
        for (int iele = 0; iele < fluiddis->NumMyRowElements(); ++iele)
        {
          DRT::Element* ele = fluiddis->lRowElement(iele);
          if (flamefront_->InterfaceHandle()->ElementCutStatus(ele->Id()) != COMBUST::InterfaceHandleCombust::uncut)
          {
            const int* nodes = ele->NodeIds();
            for (int i = 0; i < ele->NumNode(); ++i)
            {
              int gid = nodes[i];
              weights.ReplaceGlobalValues(1,&enrichednodeweight,&gid);
            }
          }
        }
      }

      // ----------------------------------------
      // loop masternodes to adjust weights of slavenodes
      // they need a small weight since they do not contribute any dofs
      // to the linear system
      {
        map<int,vector<int> >::iterator masterslavepair;

        for(masterslavepair =allcoupledcolnodes->begin();
            masterslavepair!=allcoupledcolnodes->end()  ;
            ++masterslavepair)
        {
          // get masternode
          DRT::Node*  master = fluiddis->gNode(masterslavepair->first);

          if (master->Owner() != Comm().MyPID())
            continue;

          // loop slavenodes associated with master
          for(vector<int>::iterator iter=masterslavepair->second.begin();
              iter!=masterslavepair->second.end();++iter)
          {
            int gid       =*iter;
            weights.ReplaceGlobalValues(1,&slavenodeweight,&gid);
          }
        }
      }


      // allocate graph
      RefCountPtr<Epetra_CrsGraph> nodegraph = rcp(new Epetra_CrsGraph(Copy,*noderowmap,108,false));

      // -------------------------------------------------------------
      // iterate all elements on this proc including ghosted ones
      // compute connectivity

      // standard part without master<->slave coupling
      // Note:
      // if a proc stores the appropiate ghosted elements, the resulting
      // graph will be the correct and complete graph of the distributed
      // discretization even if nodes are not ghosted.

      for (int nele = 0; nele < fluiddis->NumMyColElements(); ++nele)
      {
        // get the element
        DRT::Element* ele = fluiddis->lColElement(nele);

        // get its nodes and nodeids
        const int  nnode   = ele->NumNode();
        const int* nodeids = ele->NodeIds();

        for (int row = 0; row < nnode; ++row)
        {
          const int rownode = nodeids[row];

          // insert into line of graph only when this proc owns the node
          if (!noderowmap->MyGID(rownode))
            continue;

          // insert all neighbors from element in the graph
          for (int col = 0; col < nnode; ++col)
          {
            int colnode = nodeids[col];
            int err = nodegraph->InsertGlobalIndices(rownode,1,&colnode);
            if (err<0)
              dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
          }
        }
      }

      // -------------------------------------------------------------
      // additional coupling between master and slave
      // we do not only connect master and slave nodes but if a master/slave
      // is connected to a master/slave, we connect the corresponding slaves/master
      // as well

      for (int nele = 0; nele < fluiddis->NumMyColElements();++nele)
      {
        // get the element
        DRT::Element* ele = fluiddis->lColElement(nele);

        // get its nodes and nodeids
        const int  nnode   = ele->NumNode();
        const int* nodeids = ele->NodeIds();

        for (int row = 0; row < nnode; ++row)
        {
          const int rownode = nodeids[row];

          // insert into line of graph only when this proc owns the node
          if (!noderowmap->MyGID(rownode))
            continue;

          map<int,vector<int> >::iterator masterslavepair = allcoupledcolnodes->find(rownode);
          if(masterslavepair != allcoupledcolnodes->end())
          {
            // get all masternodes of this element
            for (int col = 0; col < nnode; ++col)
            {
              int colnode = nodeids[col];

              map<int,vector<int> >::iterator othermasterslavepair = allcoupledcolnodes->find(colnode);
              if(othermasterslavepair != allcoupledcolnodes->end())
              {
                // add connection to all slaves

                for(vector<int>::iterator iter = othermasterslavepair->second.begin();
                    iter != othermasterslavepair->second.end(); ++iter)
                {
                  int othermastersslaveindex = *iter;
                  int masterindex            = rownode;
                  int err = nodegraph->InsertGlobalIndices(rownode,1,&othermastersslaveindex);
                  if (err<0)
                    dserror("nodegraph->InsertGlobalIndices returned err=%d",err);

                  if (noderowmap->MyGID(*iter))
                  {
                    err = nodegraph->InsertGlobalIndices(*iter,1,&masterindex);
                    if (err<0)
                      dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
                  }
                }
              }
            }
          }
        }
      }

      // finalize construction of initial graph
      int err = nodegraph->FillComplete();
      if (err)
        dserror("graph->FillComplete returned %d",err);


      //--------------------------------------------------------------------------------------
      // prepare and call METIS
      //--------------------------------------------------------------------------------------

      const int myrank   = nodegraph->Comm().MyPID();
      const int numproc  = nodegraph->Comm().NumProc();

      // proc that will do the serial partitioning
      // the graph is collapsed to this proc
      // Normally this would be proc 0 but 0 always has so much to do.... ;-)
      int workrank = 1;

      // get rowmap of the graph
      const Epetra_BlockMap& tmp = nodegraph->RowMap();
      Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                        tmp.MyGlobalElements(),0,nodegraph->Comm());

      // -------------------------------------------------------------
      // build a target map that stores everything on proc workrank
      // We have arbirtary gids here and we do not tell metis about
      // them. So we have to keep rowrecv until the redistributed map is
      // build.

      // rowrecv is a fully redundant vector (size of number of nodes)
      vector<int> rowrecv(rowmap.NumGlobalElements());

      // after AllreduceEMap rowrecv contains
      //
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      // * | | .... | | * | | .... | | * ..........  * | | .... | | *
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      //   gids stored     gids stored                  gids stored
      //  on first proc  on second proc                 on last proc
      //
      // the ordering of the gids on the procs is arbitrary (as copied
      // from the map)
      LINALG::AllreduceEMap(rowrecv, rowmap);

      // construct an epetra map from the list of gids
      Epetra_Map tmap(rowmap.NumGlobalElements(),
                      // if ..........    then ............... else
                      (myrank == workrank) ? (int)rowrecv.size() : 0,
                      &rowrecv[0],
                      0,
                      rowmap.Comm());

      // export the graph to tmap
      Epetra_CrsGraph tgraph(Copy,tmap,108,false);
      Epetra_Export exporter(rowmap,tmap);
      {
        int err = tgraph.Export(*nodegraph,exporter,Add);
        if (err < 0)
          dserror("Graph export returned err=%d",err);
      }
      tgraph.FillComplete();
      tgraph.OptimizeStorage();

      // export the weights to tmap
      Epetra_Vector tweights(tmap,false);
      err = tweights.Export(weights,exporter,Insert);
      if (err < 0)
        dserror("Vector export returned err=%d",err);

      // metis requests indexes. So we need a reverse lookup from gids
      // to indexes.
      map<int,int> idxmap;
      // xadj points from index i to the index of the
      // first adjacent node
      vector<int> xadj  (rowmap.NumGlobalElements()+1);
      // a list of adjacent nodes, adressed using xadj
      vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

      // This is a vector of size n that upon successful completion stores the partition vector of the graph
      vector<int> part(tmap.NumMyElements());

      // construct reverse lookup for all procs
      for (size_t i = 0; i < rowrecv.size(); ++i)
      {
        idxmap[rowrecv[i]] = i;
      }

      if (myrank == workrank)
      {
        // ----------------------------------------

        // rowrecv(i)       rowrecv(i+1)                      node gids
        //     ^                 ^
        //     |                 |
        //     | idxmap          | idxmap
        //     |                 |
        //     v                 v
        //     i                i+1                       equivalent indices
        //     |                 |
        //     | xadj            | xadj
        //     |                 |
        //     v                 v
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //    | | | | | | | | | | | ............... | | |      adjncy
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //
        //    |       i's       |    (i+1)'s
        //    |    neighbours   |   neighbours           (numbered by equivalent indices)
        //

        int count=0;
        xadj[0] = 0;
        for (int row=0; row<tgraph.NumMyRows(); ++row)
        {
          int grid = tgraph.RowMap().GID(row);
          int numindices;
          int* lindices;
          int err = tgraph.ExtractMyRowView(row,numindices,lindices);
          if (err)
            dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);

          for (int col = 0; col < numindices; ++col)
          {
            int gcid = tgraph.ColMap().GID(lindices[col]);
            if (gcid == grid)
              continue;
            adjncy[count] = idxmap[gcid];
            ++count;
          }
          xadj[row+1] = count;
        }
      } // if (myrank == workrank)

      // broadcast xadj
      tmap.Comm().Broadcast(&xadj[0],xadj.size(),workrank);

      // broadcast adjacence (required for edge weights)
      int adjncysize = (int)adjncy.size();
      tmap.Comm().Broadcast(&adjncysize,1,workrank);
      adjncy.resize(adjncysize);
      tmap.Comm().Broadcast(&adjncy[0],adjncysize,workrank);

      // -------------------------------------------------------------
      // set a fully redundant vector of weights for edges
      vector<int> ladjwgt(adjncy.size(),0);
      vector<int>  adjwgt(adjncy.size(),0);

      for(vector<int>::iterator iter =ladjwgt.begin();
          iter!=ladjwgt.end();
          ++iter)
      {
        *iter=0;
      }

      // loop all master nodes on this proc
      {
        map<int,vector<int> >::iterator masterslavepair;

        for(masterslavepair =allcoupledcolnodes->begin();
            masterslavepair!=allcoupledcolnodes->end()  ;
            ++masterslavepair)
        {
          // get masternode
          DRT::Node*  master = fluiddis->gNode(masterslavepair->first);

          if(master->Owner()!=myrank)
          {
            continue;
          }

          map<int,int>::iterator paul = idxmap.find(master->Id());
          if (paul == idxmap.end())
            dserror("master not in reverse lookup");

          // inverse lookup
          int masterindex = idxmap[master->Id()];

          // loop slavenodes
          for(vector<int>::iterator iter = masterslavepair->second.begin();
              iter != masterslavepair->second.end(); ++iter)
          {
            DRT::Node*  slave = fluiddis->gNode(*iter);

            if(slave->Owner() != myrank)
              dserror("own master but not slave\n");

            int slaveindex = idxmap[slave->Id()];

            map<int,int>::iterator foo = idxmap.find(slave->Id());
            if (foo == idxmap.end())
              dserror("slave not in reverse lookup");

            // -------------------------------------------------------------
            // connections between master and slavenodes are very strong
            // we do not want to partition between master and slave nodes
            for(int j = xadj[masterindex]; j < xadj[masterindex+1]; ++j)
            {
              if(adjncy[j] == slaveindex)
              {
                ladjwgt[j] = 100;
              }
            }

            for(int j = xadj[slaveindex]; j < xadj[slaveindex+1]; ++j)
            {
              if(adjncy[j] == masterindex)
              {
                ladjwgt[j] = 100;
              }
            }
          }
        }
      }

      // do communication to aquire edge weight information from all procs
      tmap.Comm().SumAll(&ladjwgt[0], &adjwgt[0], adjwgt.size());

      // the standard edge weight is one
      for(vector<int>::iterator iter =adjwgt.begin();
          iter!=adjwgt.end();
          ++iter)
      {
        if(*iter==0)
        *iter=1;
      }

      // the reverse lookup is not required anymore
      idxmap.clear();

      // -------------------------------------------------------------
      // do partitioning using metis on workrank
      if (myrank == workrank)
      {
        // the vertex weights
        vector<int> vwgt(tweights.MyLength());
        for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

        // 0 No weights (vwgts and adjwgt are NULL)
        // 1 Weights on the edges only (vwgts = NULL)
        // 2 Weights on the vertices only (adjwgt = NULL)
        // 3 Weights both on vertices and edges.
        int wgtflag=3;
        // 0 C-style numbering is assumed that starts from 0
        // 1 Fortran-style numbering is assumed that starts from 1
        int numflag=0;
        // The number of parts to partition the graph.
        int npart=numproc;
        // This is an array of 5 integers that is used to pass parameters for the various phases of the algorithm.
        // If options[0]=0 then default values are used. If options[0]=1, then the remaining four elements of
        // options are interpreted as follows:
        // options[1]    Determines matching type. Possible values are:
        //               1 Random Matching (RM)
        //               2 Heavy-Edge Matching (HEM)
        //               3 Sorted Heavy-Edge Matching (SHEM) (Default)
        //               Experiments has shown that both HEM and SHEM perform quite well.
        // options[2]    Determines the algorithm used during initial partitioning. Possible values are:
        //               1 Region Growing (Default)
        // options[3]    Determines the algorithm used for re%Gï¬%@nement. Possible values are:
        //               1 Early-Exit Boundary FM re%Gï¬%@nement (Default)
        // options[4]    Used for debugging purposes. Always set it to 0 (Default).
        int options[5] = { 0,3,1,1,0 };
        // Upon successful completion, this variable stores the number of edges that are cut by the partition.
        int edgecut=0;
        // The number of vertices in the graph.
        int nummyele = tmap.NumMyElements();

        cout << "proc " <<  myrank << " repartition graph using metis" << endl;

        if (numproc<8) // better for smaller no. of partitions
        {
          METIS_PartGraphRecursive(&nummyele,
                                   &xadj[0],
                                   &adjncy[0],
                                   &vwgt[0],
                                   &adjwgt[0],
                                   &wgtflag,
                                   &numflag,
                                   &npart,
                                   options,
                                   &edgecut,
                                   &part[0]);

          cout << "METIS_PartGraphRecursive produced edgecut of " << edgecut << endl;
        }
        else
        {
          METIS_PartGraphKway(&nummyele,
                              &xadj[0],
                              &adjncy[0],
                              &vwgt[0],
                              &adjwgt[0],
                              &wgtflag,
                              &numflag,
                              &npart,
                              options,
                              &edgecut,
                              &part[0]);
          cout << "METIS_PartGraphKway produced edgecut of " << edgecut << endl;
        }
      } // if (myrank==workrank)

      // broadcast partitioning result
      int size = tmap.NumMyElements();
      tmap.Comm().Broadcast(&size,1,workrank);
      part.resize(size);
      tmap.Comm().Broadcast(&part[0],size,workrank);

      // loop part and count no. of nodes belonging to me
      // (we reuse part to save on memory)
      int count = 0;
      for (int i = 0; i < size; ++i)
      {
        if (part[i]==myrank)
        {
          part[count] = rowrecv[i];
          ++count;
        }
      }

      // rowrecv is done
      rowrecv.clear();

      // create map with new layout
      Epetra_Map newmap(size,count,&part[0],0,nodegraph->Comm());

      // create the new graph and export to it
      RefCountPtr<Epetra_CrsGraph> newnodegraph;

      newnodegraph = rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
      Epetra_Export exporter2(nodegraph->RowMap(),newmap);
      err = newnodegraph->Export(*nodegraph,exporter2,Add);
      if (err<0)
        dserror("Graph export returned err=%d",err);
      newnodegraph->FillComplete();
      newnodegraph->OptimizeStorage();


      //--------------------------------------------------------------------------------------
      // Ensure the new distribution is valid
      //--------------------------------------------------------------------------------------
      {
        // the rowmap will become the new distribution of nodes
        const Epetra_BlockMap rntmp = newnodegraph->RowMap();
        Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,Comm());

        // the column map will become the new ghosted distribution of nodes
        const Epetra_BlockMap Mcntmp = newnodegraph->ColMap();
        Epetra_Map newnodecolmap(-1,Mcntmp.NumMyElements(),Mcntmp.MyGlobalElements(),0,Comm());

        Teuchos::RCP<Epetra_Map> elerowmap;
        Teuchos::RCP<Epetra_Map> elecolmap;

        ScaTraField().Discretization()->BuildElementRowColumn(newnoderowmap,newnodecolmap,elerowmap,elecolmap);

        int myrowele = elerowmap->NumMyElements();
        int minrowele = 1;
        Comm().MinAll(&myrowele, &minrowele, 1);
        if (minrowele <= 0)
        {
          if (Comm().MyPID() == 0)
            std::cout << "A processor would be left without row elements --> Redistribution aborted\n"
                      << "----------------------------------------------------------------------------" << std::endl;
          return;
        }
      }


      //--------------------------------------------------------------------------------------
      // call ScaTra and Fluid field and pass the new graph, so they can manage all
      // redistribution themselves
      //--------------------------------------------------------------------------------------
      if(Comm().MyPID()==0)
        cout << "Redistributing ScaTra Discretization                                ... " << flush;

      ScaTraField().Redistribute(newnodegraph);

      if (reinitaction_ == INPAR::COMBUST::reinitaction_sussman)
      {
        if(Comm().MyPID()==0)
          cout << "done\nRedistributing ScaTra Reinit Discretization                         ... " << flush;
        dserror("Implement a reinit_pde_->Redistribute() method for pde_reinitialization class if necessary");
      }

      if(Comm().MyPID()==0)
        cout << "Redistributing Fluid Discretization                                 ... " << flush;

      FluidField().Redistribute(newnodegraph);

      if(Comm().MyPID()==0)
        cout << "done\nUpdating interface                                                  ... " << flush;

      // update flame front according to evolved G-function field
      // remark: for only one FGI iteration, 'phinpip_' == ScaTraField().Phin()
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField().Phin(), ScaTraField().Phinp());

      if(Comm().MyPID()==0)
        cout << "Transfering state vectors to new distribution                       ... " << flush;

      FluidField().TransferVectorsToNewDistribution(flamefront_);

      if(Comm().MyPID()==0)
        cout << "done" << endl;

      //--------------------------------------------------------------------------------------
      // with the scatra- and fluid-field updated it is time to do the same with the algorithm
      //--------------------------------------------------------------------------------------
      Teuchos::RCP<Epetra_Vector> old;

      if (velnpi_ != Teuchos::null)
      {
        old = velnpi_;
        velnpi_ = rcp(new Epetra_Vector(FluidField().StdVelnp()->Map()),true);
        LINALG::Export(*old, *velnpi_);
      }

      if (phinpi_ != Teuchos::null)
      {
        old = phinpi_;
        phinpi_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
        LINALG::Export(*old, *phinpi_);
      }

      //-----------------------------------------------------------------------------
      // redistibute vectors holding the means in the general mean statistics manager
      //-----------------------------------------------------------------------------
      // if applicable, provide scatra data to the turbulence statistics
      if (FluidField().TurbulenceStatisticManager() != Teuchos::null)
      {
//        if(Comm().MyPID()==0)
//          cout << "Updating pointer to ScaTra vector                                               ... " << flush;

        // update the statistics manager to the new ScaTra discretization
//        FluidField().TurbulenceStatisticManager()
//            ->AddScaTraResults(ScaTraField().Discretization(),ScaTraField().Phinp());

        if(Comm().MyPID()==0)
          cout << "Redistributing General Mean Statistics Manager                      ... " << flush;

//        // redistribute with redistributed standard fluid dofset
        if ( FluidField().TurbulenceStatisticManager()->GetTurbulenceStatisticsGeneralMean()!=Teuchos::null)
          FluidField().TurbulenceStatisticManager()->GetTurbulenceStatisticsGeneralMean()
              ->Redistribute(FluidField().DofSet(), FluidField().Discretization());

        if(Comm().MyPID()==0)
          cout << "done" << endl;
      }

      if (Comm().MyPID() == 0)
        cout << "----------------------------------------------------------------------------" << endl;

      // make sure the redistribution ratio will be reset
      evaltimeratio_ = 1e12;
    }
  }
  return;
}
