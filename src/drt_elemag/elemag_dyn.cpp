/*----------------------------------------------------------------------*/
/*! \file

\brief Main control routine for electromagnetic simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "elemag_dyn.H"
#include "elemag_timeint.H"
#include "elemag_ele.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_elemag.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
//#include "../linalg/linalg_utils_sparse_algebra_print.H"

// Clone discreization
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_elemag/elemag_utils_clonestrategy.H"

// Scatra init
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_stat_hdg.H"

void redistributescatradis(
    Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<DRT::Discretization> elemagdis)
{
  const Teuchos::ParameterList& elemagparams = DRT::Problem::Instance()->ElectromagneticParams();
  bool meshconform = DRT::INPUT::IntegralValue<bool>(elemagparams, "MESHCONFORM");
  if (meshconform == false) dserror("at the moment only implemented for conforming meshes.");

  // redistribute scatra discretization, since it usually has way less dofs than acou and does not
  // need as many processors
  // if(DRT::INPUT::IntegralValue<bool>(elemagparams,"REDISTRIBUTESCATRA"))
  {
    std::vector<int> rownodeids;
    std::vector<int> colnodeids;

    int minelemagnodegid = elemagdis->NodeRowMap()->MinAllGID();
    int minscanodegid = scatradis->NodeRowMap()->MinAllGID();

    for (int i = minscanodegid; i < minscanodegid + scatradis->NumGlobalNodes(); ++i)
    {
      // check who owns the eelctromagnetic node
      int lidrownode = elemagdis->NodeRowMap()->LID(i - minscanodegid + minelemagnodegid);
      int lidcolnode = elemagdis->NodeColMap()->LID(i - minscanodegid + minelemagnodegid);

      if (lidrownode >= 0) rownodeids.push_back(i);
      if (lidcolnode >= 0) colnodeids.push_back(i);
    }

    // create new maps
    Teuchos::RCP<Epetra_Map> newnoderowmap = Teuchos::rcp(new Epetra_Map(
        scatradis->NumGlobalNodes(), rownodeids.size(), &rownodeids[0], 0, elemagdis->Comm()));
    Teuchos::RCP<Epetra_Map> newnodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, colnodeids.size(), &colnodeids[0], 0, elemagdis->Comm()));

    // redistribute
    scatradis->Redistribute(*newnoderowmap, *newnodecolmap, true, true, true, true, true);
    std::cout << "warning: this function is not tested for empty processors or processors which "
                 "accidently only get one or two elements, test it and implement something better"
              << std::endl;
  }

  return;
}

void electromagnetics_drt()
{
  // declare abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // The function NumDofPerElementAuxiliary() of the electromagnetic elements return nsd_*2. This
  // does not assure that the code will work in any case (more spatial dimensions might give
  // problems)
  if (problem->NDim() != 3)
  {
    dserror(
        "The implementation of electromagnetic propagation only supports 3D problems.\n"
        "It is necessary to change the spatial dimension of your problem.");
  }

  // declare problem-specific parameter list for electromagnetics
  const Teuchos::ParameterList& elemagparams = problem->ElectromagneticParams();

  // declare discretization and check their existence
  Teuchos::RCP<DRT::DiscretizationHDG> elemagdishdg =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationHDG>(problem->GetDis("elemag"));
  if (elemagdishdg == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

#ifdef DEBUG
  elemagdishdg->PrintFaces(std::cout);
#endif

  // declare communicator and print module information to screen
  const Epetra_Comm& comm = elemagdishdg->Comm();
  if (comm.MyPID() == 0)
  {
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "---------- You are now about to enter the module for electromagnetics! ----------"
              << std::endl;
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
  }

  // call fill complete on discretization
  if (not elemagdishdg->Filled() || not elemagdishdg->HaveDofs()) elemagdishdg->FillComplete();
  // Asking the discretization how many internal DOF the elements have and creating the additional
  // DofSet
  int eledofs = dynamic_cast<DRT::ELEMENTS::Elemag*>(elemagdishdg->lColElement(0))
                    ->NumDofPerElementAuxiliary();
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0, eledofs, 0, false));
  elemagdishdg->AddDofSet(dofsetaux);

  // call fill complete on discretization
  elemagdishdg->FillComplete();

  // create solver
  const int linsolvernumber_elemag = elemagparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber_elemag == (-1))
    dserror(
        "There is not any linear solver defined for electromagnetic problem. Please set "
        "LINEAR_SOLVER in ELECTROMAGNETIC DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      problem->SolverParams(linsolvernumber_elemag), comm, problem->ErrorFile()->Handle()));

  // declare output writer
  Teuchos::RCP<IO::DiscretizationWriter> output = elemagdishdg->Writer();

  // declare electromagnetic parameter list
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList(elemagparams));

  // set restart step if required
  int restart = problem->Restart();
  params->set<int>("restart", restart);

  // create algorithm depending on time-integration scheme
  INPAR::ELEMAG::DynamicType elemagdyna =
      DRT::INPUT::IntegralValue<INPAR::ELEMAG::DynamicType>(elemagparams, "TIMEINT");
  Teuchos::RCP<ELEMAG::ElemagTimeInt> elemagalgo;
  switch (elemagdyna)
  {
    case INPAR::ELEMAG::elemag_ost:
    {
      dserror("One step theta not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimIntOST(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_bdf1:
    case INPAR::ELEMAG::elemag_bdf2:
    case INPAR::ELEMAG::elemag_bdf4:
    {
      elemagalgo = Teuchos::rcp(new ELEMAG::ElemagTimeInt(elemagdishdg, solver, params, output));
      break;
    }
    case INPAR::ELEMAG::elemag_genAlpha:
    {
      dserror("Generalized-alpha method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::ElemagGenAlpha(elemagdishdg, solver, params,
      // output));
      break;
    }
    case INPAR::ELEMAG::elemag_explicit_euler:
    {
      dserror("Explicit euler method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntExplEuler(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_rk:
    {
      dserror("Runge-Kutta methods not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntRK(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_cn:
    {
      dserror("Crank-Nicolson method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntCN(elemagdishdg,solver,params,output));
      break;
    }
    default:
      dserror("Unknown time-integration scheme for problem type electromagnetics");
      break;
  }

  // Initialize the evolution algorithm
  elemagalgo->Init();

  // set initial field
  if (restart)
    elemagalgo->ReadRestart(restart);
  else
  {
    INPAR::ELEMAG::InitialField init =
        DRT::INPUT::IntegralValue<INPAR::ELEMAG::InitialField>(elemagparams, "INITIALFIELD");

    bool ishdg = false;

    switch (init)
    {
      case INPAR::ELEMAG::initfield_scatra_hdg:
        ishdg = true;
      case INPAR::ELEMAG::initfield_scatra:
      {
        Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(elemagdishdg->Comm().Clone());

        Teuchos::RCP<DRT::Discretization> scatradis;

        if (ishdg)
        {
          scatradis = Teuchos::rcp(new DRT::DiscretizationHDG((std::string) "scatra", newcomm));

          scatradis->FillComplete();

          DRT::UTILS::CloneDiscretization<
              ELEMAG::UTILS::ScatraCloneStrategy<ShapeFunctionType::shapefunction_hdg>>(
              elemagdishdg, scatradis);
        }
        else
        {
          scatradis = Teuchos::rcp(new DRT::Discretization((std::string) "scatra", newcomm));
          scatradis->FillComplete();

          DRT::UTILS::CloneDiscretization<
              ELEMAG::UTILS::ScatraCloneStrategy<ShapeFunctionType::shapefunction_polynomial>>(
              elemagdishdg, scatradis);
        }

        // call fill complete on discretization
        scatradis->FillComplete();

        Teuchos::RCP<IO::DiscretizationWriter> output_scatra = scatradis->Writer();

        if (0)
        {
          redistributescatradis(scatradis, elemagdishdg);
          scatradis->FillComplete();
        }

        // This is necessary to have the dirichlet conditions done also in the scatra problmem. It
        // might be necessary to rethink how things are handled inside the
        // DRT::UTILS::DbcHDG::DoDirichletCondition.
        problem->SetProblemType(prb_scatra);

        // access the problem-specific parameter list
        const Teuchos::ParameterList& scatradyn =
            DRT::Problem::Instance()->ScalarTransportDynamicParams();

        // do the scatra
        const INPAR::SCATRA::VelocityField veltype =
            DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");
        switch (veltype)
        {
          case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
          {
            // we directly use the elements from the scalar transport elements section
            if (scatradis->NumGlobalNodes() == 0)
              dserror("No elements in the ---TRANSPORT ELEMENTS section");

            // add proxy of velocity related degrees of freedom to scatra discretization
            Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
                Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(
                    DRT::Problem::Instance()->NDim() + 1, 0, 0, true));
            if (scatradis->AddDofSet(dofsetaux) != 1)
              dserror("Scatra discretization has illegal number of dofsets!");

            // finalize discretization
            scatradis->FillComplete(true, true, true);

            // create scatra output
            // Teuchos::RCP<IO::DiscretizationWriter> scatraoutput = scatradis->Writer();
            // output_scatra->NewResultFile("scatra_file");
            // output->WriteMesh(0, 0.0);
            // access the problem-specific parameter list
            Teuchos::RCP<Teuchos::ParameterList> scatraparams =
                Teuchos::rcp(new Teuchos::ParameterList(
                    DRT::Problem::Instance()->ScalarTransportDynamicParams()));

            // TODO (berardocco) Might want to add the scatra section in the input file to avoid
            // adding params to the elemag or using existing ones for scatra purposes
            scatraparams->set("TIMEINTEGR", "Stationary");
            scatraparams->set("NUMSTEP", 1);
            // This way we avoid writing results and restart
            scatraparams->set("RESULTSEVRY", 1000);
            scatraparams->set("RESTARTEVRY", 1000);
            // This has to be changed accordingly to the intial time
            // As of now elemag simulation can only start at 0.
            scatraparams->set<double>("MAXTIME", 0 + elemagparams.get<double>("TIMESTEP"));
            scatraparams->set("TIMESTEP", elemagparams.get<double>("TIMESTEP"));
            // scatraparams->set("CALCERROR", "error_by_function");
            // scatraparams->set("CALCERRORNO", elemagparams.get<int>("ERRORFUNCNO"));

            // The solver type still has to be fixed as the problem is linear but the steady state
            // does not always behave correctly with linear solvers.
            scatraparams->set("SOLVERTYPE", "nonlinear");

            // create necessary extra parameter list for scatra
            //{
            Teuchos::RCP<Teuchos::ParameterList> scatraextraparams;
            scatraextraparams = Teuchos::rcp(new Teuchos::ParameterList());
            scatraextraparams->set<FILE*>(
                "err file", DRT::Problem::Instance()->ErrorFile()->Handle());
            scatraextraparams->set<bool>("isale", false);
            const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
            scatraextraparams->sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");
            scatraextraparams->sublist("SUBGRID VISCOSITY") = fdyn.sublist("SUBGRID VISCOSITY");
            scatraextraparams->sublist("MULTIFRACTAL SUBGRID SCALES") =
                fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
            scatraextraparams->sublist("TURBULENT INFLOW") = fdyn.sublist("TURBULENT INFLOW");
            //}

            scatraextraparams->set("ELECTROMAGNETICDIFFUSION", true);
            scatraextraparams->set("EMDSOURCE", elemagparams.get<int>("SOURCEFUNCNO"));

            // In case the scatra solver is not defined just use the elemag one
            if (scatraparams->get<int>("LINEAR_SOLVER") == -1)
              scatraparams->set<int>("LINEAR_SOLVER", elemagparams.get<int>("LINEAR_SOLVER"));

            // create solver
            Teuchos::RCP<LINALG::Solver> scatrasolver = Teuchos::rcp(new LINALG::Solver(
                DRT::Problem::Instance()->SolverParams(scatraparams->get<int>("LINEAR_SOLVER")),
                scatradis->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

            // create instance of scalar transport basis algorithm (empty fluid discretization)
            Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatraalgo;
            if (ishdg)
            {
              // Add parameters for HDG
              scatraparams->sublist("STABILIZATION").set("STABTYPE", "centered");
              scatraparams->sublist("STABILIZATION").set("DEFINITION_TAU", "Numerical_Value");
              // If the input file does not specify a tau parameter then use the one given to the
              // elemag discretization
              if (scatraparams->sublist("STABILIZATION").get<double>("TAU_VALUE") == 0.0)
                scatraparams->sublist("STABILIZATION")
                    .set("TAU_VALUE", elemagparams.get<double>("TAU"));

              scatraalgo = Teuchos::rcp(new SCATRA::TimIntStationaryHDG(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }
            else
            {
              // Add parameters for CG
              // There is no need for stabilization as the problem is a pure diffusion problem
              scatraparams->sublist("STABILIZATION").set("STABTYPE", "no_stabilization");
              scatraalgo = Teuchos::rcp(new SCATRA::TimIntStationary(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }

            // scatraparams->print(std::cout);

            scatraalgo->Init();
            scatraalgo->Setup();
            scatraalgo->SetVelocityField(1);
            scatraalgo->TimeLoop();

            // scatraalgo->ComputeInteriorValues();

            // Create a vector that is going to be filled differently depending on the
            // discretization. If HDG we already have the information about the gradient, otherwise
            // the gradient has to be computed.
            Teuchos::RCP<Epetra_Vector> phi;

            // If it is an HDG discretization return the interior variables else return the nodal
            // values
            if (ishdg)
              phi = Teuchos::rcp_dynamic_cast<SCATRA::TimIntStationaryHDG>(scatraalgo)
                        ->ReturnIntPhinp();
            else
              phi = scatraalgo->Phinp();

            // This is a shortcut for output reason
            // TODO (berardocco) Fix the output
            output->CreateNewResultAndMeshFile();

            // Given the results of the scatra solver obtain the initial value of the electric field
            elemagalgo->SetInitialElectricField(phi, scatradis);

            // Once work is done change back to problem elemag
            problem->SetProblemType(prb_elemag);

            break;
          }
          default:
            dserror(
                "Does not make sense to have a velocity field to initialize the electric potential "
                "field.\nCheck your input file.");
            break;
        }
        break;
      }
      default:
      {
        int startfuncno = elemagparams.get<int>("STARTFUNCNO");
        elemagalgo->SetInitialField(init, startfuncno);
        break;
      }
    }
  }

  // print information to screen
  elemagalgo->PrintInformationToScreen();

  elemagalgo->Integrate();

  // PrintMatrixInMatlabFormat("matrix_matlab");
  // LINALG::PrintMatrixInMatlabFormat("hdg_elemag", elemagalgo->SystemMatrix());
  // elemagalgo->SpySysmat("matrix.mat");

  // Computing the error at the las time step (the conditional stateme nt is inside for now)
  if (DRT::INPUT::IntegralValue<bool>(elemagparams, "CALCERR"))
  {
    Teuchos::RCP<Epetra_SerialDenseVector> errors = elemagalgo->ComputeError();
    elemagalgo->PrintErrors(errors);
  }

  // print computing time
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // do result test if required
  problem->AddFieldTest(elemagalgo->CreateFieldTest());
  problem->TestAll(comm);

  return;
}
