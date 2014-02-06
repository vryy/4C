/*!----------------------------------------------------------------------
\file acou_dyn.cpp
\brief Main control routine for acoustic simulations

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "acou_dyn.H"
#include "acou_impl.H"
#include "acou_impl_euler.H"
#include "acou_impl_trap.H"
#include "acou_impl_bdf.H"
#include "acou_impl_dirk.H"
#include "acou_inv_analysis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_acou.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_comm/comm_utils.H"

void printacoulogo()
{
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "-----------------   Welcome to the problem type Acoustics       -----------------" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  return;
}

void printacouinvlogo()
{
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "-----------------   Welcome to the problem type Acoustics       -----------------" << std::endl;
  std::cout << "-----------------            Inverse analysis                   -----------------" << std::endl;
  std::cout << "-----------------    Photoacoustic image reconstruction         -----------------" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  return;
}

void acoustics_drt()
{
  // get input lists
  const Teuchos::ParameterList& acouparams = DRT::Problem::Instance()->AcousticParams();

  //if(DRT::Problem::Instance()->SpatialApproximation() != "HDG") dserror("you have to set SHAPEFCT in parameter list PROBLEM TYP to HDG for Acoustics");

  bool invanalysis = (DRT::INPUT::IntegralValue<INPAR::ACOU::InvAnalysisType>(acouparams,"INV_ANALYSIS") != INPAR::ACOU::inv_none);

  // access the discretization
  Teuchos::RCP<DRT::DiscretizationHDG> acoudishdg = Teuchos::rcp_static_cast<DRT::DiscretizationHDG>(DRT::Problem::Instance()->GetDis("acou"));
  if (acoudishdg == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

  const int dim = DRT::Problem::Instance()->NDim();
  int nscalardofs = 1;
  for (int i=0; i<dim; ++i)
    nscalardofs *= 4;
  const int elementndof =  dim * nscalardofs       // velocity DoFs
                           +
                           nscalardofs;            // pressure DoFs
  // set degrees of freedom in the discretization
//  Teuchos::RCP<DRT::IndependentDofSet> secondary = Teuchos::rcp(new DRT::IndependentDofSet());
//  acoudishdg->AddDofSet(secondary);
  acoudishdg->BuildDofSetAuxProxy(0,elementndof,false);

  // call fill complete on acoustical discretization
  if (not acoudishdg->Filled() || not acoudishdg->HaveDofs()) acoudishdg->FillComplete();

  // print problem specific logo
  if(!acoudishdg->Comm().MyPID()) { if (!invanalysis) printacoulogo(); else printacouinvlogo(); }

  const int linsolvernumber_acou = acouparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber_acou == (-1))
    dserror("no linear solver defined for acoustical problem. Please set LINEAR_SOLVER in ACOUSTIC DYNAMIC to a valid number!");

  // create solver
  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber_acou),
                                    acoudishdg->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
  if ( DRT::Problem::Instance()->SolverParams(linsolvernumber_acou).get<int>("AZREUSE") < 10 )
    dserror("note: you can set AZREUSE to NUMSTEP (ACOUSTIC DYNAMIC) because the system matrix does not change -> save lots of CPU time");

  // create output
  Teuchos::RCP<IO::DiscretizationWriter> output = acoudishdg->Writer();

  // create acoustical parameter list
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList(acouparams));
  if(invanalysis)
  {
    params->set<bool>("invana",true);
  }
  else
  {
    params->set<bool>("invana",false);
    params->set<bool>("adjoint",false);
  }

  // set pulse duration
  double pulse = acouparams.get<double>("PULSEDURATION");

  // set restart step if we do not perform inverse analysis
  int restart = -1;
  if(!invanalysis)
  {
    restart = DRT::Problem::Instance()->Restart();
  }
  params->set<int>("restart",restart);


  // the main part of this function:
  if(!invanalysis)
  {
    // do we do a PHOTOacoustic simulation? if so -> do a SCATRA simulation for initial pressure distribution
    bool photoacou = DRT::INPUT::IntegralValue<bool>(acouparams,"PHOTOACOU");

    // create algorithm depending on the temporal discretization
    Teuchos::RCP<ACOU::AcouImplicitTimeInt> acoualgo;

    INPAR::ACOU::DynamicType dyna = DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(acouparams,"TIMEINT");
    switch(dyna)
    {
    case INPAR::ACOU::acou_impleuler:
    {
      acoualgo = Teuchos::rcp(new ACOU::TimIntImplEuler(acoudishdg,solver,params,output));
      break;
    }
    case INPAR::ACOU::acou_trapezoidal:
    {
      acoualgo = Teuchos::rcp(new ACOU::TimIntImplTrap(acoudishdg,solver,params,output));
      break;
    }
    case INPAR::ACOU::acou_bdf2:
    case INPAR::ACOU::acou_bdf3:
    case INPAR::ACOU::acou_bdf4:
    {
      acoualgo = Teuchos::rcp(new ACOU::TimIntImplBDF(acoudishdg,solver,params,output));
      break;
    }
    case INPAR::ACOU::acou_dirk23:
    case INPAR::ACOU::acou_dirk33:
    case INPAR::ACOU::acou_dirk34:
    case INPAR::ACOU::acou_dirk54:
    {
      acoualgo = Teuchos::rcp(new ACOU::TimIntImplDIRK(acoudishdg,solver,params,output));
      break;
    }
    default:
      dserror("Unknown time integration scheme for problem type Acoustics");
      break;
    }

    // in case we do a photoacoustic simulation
    Teuchos::RCP<Epetra_Vector> scatrainitialpress;

    if ( photoacou  && restart == 0)
    {
      // access the problem-specific parameter list
      const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

      // access the scatra discretization
      Teuchos::RCP<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");
      scatradis->FillComplete();

      if(scatradis->NumGlobalElements()==0)
        dserror("you said you want to do photoacoustics but you did not supply TRANSP elements");

      // set velocity field
      const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");
      switch (veltype)
      {
        case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
        case INPAR::SCATRA::velocity_function:  // function
        {
          // we directly use the elements from the scalar transport elements section
          if (scatradis->NumGlobalNodes()==0)
            dserror("No elements in the ---TRANSPORT ELEMENTS section");

          // get linear solver id from SCALAR TRANSPORT DYNAMIC
          const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
          if (linsolvernumber == (-1))
            dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

          // create instance of scalar transport basis algorithm (empty fluid discretization)
          Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,false,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber)));

          // set velocity field
          //(this is done only once. Time-dependent velocity fields are not supported)
          (scatraonly->ScaTraField())->SetVelocityField();

          // enter time loop to solve problem with given convective velocity
          (scatraonly->ScaTraField())->TimeLoop();

          // get the solution of the scatra problem
          scatrainitialpress = (scatraonly->ScaTraField())->Phinp();

          break;
        }
        default:
          dserror("unknown velocity field type for transport of passive scalar in problem type Acoustics");
          break;
      }
    }

    // set initial field
    if (restart) // standard restart scenario
    {
      acoualgo->ReadRestart(restart);
    }
    else if (!photoacou) // standard acoustic simulation
    {
      int startfuncno = acouparams.get<int>("STARTFUNCNO");
      acoualgo->SetInitialField(startfuncno,pulse);
    }
    else // initial field calculated in dependence on light distribution - photoacoustic setting
    {
      bool meshconform = DRT::INPUT::IntegralValue<bool>(acouparams,"MESHCONFORM");
      // calculate initial pressure distribution by means of the scalar transport problem
      acoualgo->SetInitialPhotoAcousticField(pulse,scatrainitialpress,DRT::Problem::Instance()->GetDis("scatra"), meshconform);
    }
    acoualgo->Integrate();

    // print monitoring of time consumption
    Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(acoudishdg->Comm());
    Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

    DRT::Problem::Instance()->AddFieldTest(acoualgo->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(acoudishdg->Comm());
  } // if(!invanalysis)
  else
  {
    // error message in case someone does not exactly know, what to do
    if ( !DRT::INPUT::IntegralValue<bool>(acouparams,"PHOTOACOU") )
      dserror("you chose to do a photoacoustic image reconstruction, so please set PHOTOACOU to Yes");

    // access the problem-specific parameter list
    const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

    // access the scatra discretization
    Teuchos::RCP<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");

    // ensure that all dofs are assigned in the right order
    scatradis->FillComplete();

    if ( scatradis->NumGlobalElements() == 0 )
      dserror("you said you want to do photoacoustics but you did not supply TRANSP elements");

    // till here it was all the same...
    // now, we set up the inverse analysis algorithm
    const int linsolvernumber_scatra = scatradyn.get<int>("LINEAR_SOLVER");
    if ( linsolvernumber_scatra == (-1) )
      dserror("no linear solver defined for acoustical problem. Please set LINEAR_SOLVER in ACOUSTIC DYNAMIC to a valid number!");

    ACOU::InvAnalysis myinverseproblem(scatradis,acoudishdg,scatradyn,DRT::Problem::Instance()->SolverParams(linsolvernumber_scatra),params,solver,output);
    myinverseproblem.Integrate();

    // do testing if required
    DRT::Problem::Instance()->AddFieldTest(myinverseproblem.CreateFieldTest());
    DRT::Problem::Instance()->TestAll(scatradis->Comm());

  } // else ** if(!invanalysis)
  return;

}

