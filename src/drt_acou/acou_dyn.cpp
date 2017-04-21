/*!----------------------------------------------------------------------
\file acou_dyn.cpp
\brief Main control routine for acoustic simulations
\level 2
<pre>
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "acou_dyn.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"
#include "acou_timeint.H"
#include "acou_impl_euler.H"
#include "pat_imagereconstruction.H"
#include "acou_expl.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_acou.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_scatra/scatra_timint_implicit.H"

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

void redistributescatradis(Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<DRT::Discretization> acoudis)
{
  const Teuchos::ParameterList& acouparams = DRT::Problem::Instance()->AcousticParams();
  bool meshconform = DRT::INPUT::IntegralValue<bool>(acouparams,"MESHCONFORM");
  if(meshconform==false)
    dserror("at the moment only implemented for conforming meshes.");

  // redistribute scatra discretization, since it usually has way less dofs than acou and does not need as many processors
  //if(DRT::INPUT::IntegralValue<bool>(acouparams,"REDISTRIBUTESCATRA"))
  {
    std::vector<int> rownodeids;
    std::vector<int> colnodeids;

    int minacounodegid = acoudis->NodeRowMap()->MinAllGID();
    int minscanodegid = scatradis->NodeRowMap()->MinAllGID();

    for(int i=minscanodegid; i<minscanodegid+scatradis->NumGlobalNodes(); ++i)
    {
      // check who owns the acoustical node
      int lidrownode = acoudis->NodeRowMap()->LID(i-minscanodegid+minacounodegid);
      int lidcolnode = acoudis->NodeColMap()->LID(i-minscanodegid+minacounodegid);

      if(lidrownode>=0)
        rownodeids.push_back(i);
      if(lidcolnode>=0)
        colnodeids.push_back(i);
    }

    // create new maps
    Teuchos::RCP<Epetra_Map> newnoderowmap = Teuchos::rcp(new Epetra_Map(scatradis->NumGlobalNodes(),rownodeids.size(),&rownodeids[0],0,acoudis->Comm()));
    Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,colnodeids.size(),&colnodeids[0],0,acoudis->Comm()));

    // redistribute
    scatradis->Redistribute(*newnoderowmap,*newnodecolmap,true,true,true,true,true);
    std::cout<<"warning: this function is not tested for empty processors or processors which accidently only get one or two elements, test it and implement something better"<<std::endl;
  }

  return;
}

void acoustics_drt()
{
  // get input lists
  const Teuchos::ParameterList& acouparams = DRT::Problem::Instance()->AcousticParams();

  // do we want to do inverse analysis?
  bool invanalysis = (DRT::INPUT::IntegralValue<INPAR::ACOU::InvAnalysisType>(acouparams,"INV_ANALYSIS") != INPAR::ACOU::pat_none);

  // access the discretization
  Teuchos::RCP<DRT::DiscretizationHDG> acoudishdg = Teuchos::rcp_dynamic_cast<DRT::DiscretizationHDG>(DRT::Problem::Instance()->GetDis("acou"));
  if (acoudishdg == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

  // call fill complete on acoustical discretization
  if (not acoudishdg->Filled() || not acoudishdg->HaveDofs()) acoudishdg->FillComplete();

  // for explicit time integration, need to provide some additional ghosting
  {
    switch(DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(acouparams,"TIMEINT"))
    {
    case INPAR::ACOU::acou_expleuler:
    case INPAR::ACOU::acou_classrk4:
    case INPAR::ACOU::acou_lsrk45reg2:
    case INPAR::ACOU::acou_lsrk33reg2:
    case INPAR::ACOU::acou_lsrk45reg3:
    case INPAR::ACOU::acou_ssprk:
    case INPAR::ACOU::acou_ader:
    case INPAR::ACOU::acou_ader_lts:
    {
      acoudishdg->AddElementGhostLayer();
      break;
    }
    default:
      break;
    }
  }

  // if explicit time integration and float evaluation, we need to set something for denormal numbers:
  {
    if(DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(acouparams,"TIMEINT") != INPAR::ACOU::acou_impleuler)
      if( (DRT::INPUT::IntegralValue<bool>(acouparams,"DOUBLEORFLOAT")) == false)
      {
        // change mode for rounding: denormals are flushed to zero to avoid computing
        // on denormals which can slow down things.
        // TODO: uses syntax specific to the gcc compiler, must avoid that other
        // compilers end up here. ICC is particularly tricky because it exports
        // the __GNUC__ symbol.
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
        #define MXCSR_DAZ (1 << 6)      /* Enable denormals are zero mode */
        #define MXCSR_FTZ (1 << 15)     /* Enable flush to zero mode */

        unsigned int mxcsr = __builtin_ia32_stmxcsr ();
        mxcsr |= MXCSR_DAZ | MXCSR_FTZ;
        __builtin_ia32_ldmxcsr (mxcsr);
#endif
      }
  }

  // set degrees of freedom in the discretization
  //acoudishdg->BuildDofSetAuxProxy(0,elementndof,0,false);
  // build map
  const Teuchos::RCP<Epetra_IntVector> eledofs = Teuchos::rcp(new Epetra_IntVector(*acoudishdg->ElementColMap()));
  for(int i=0; i<acoudishdg->NumMyColElements(); ++i)
  {
    (*eledofs)[i] = dynamic_cast<DRT::ELEMENTS::Acou*>(acoudishdg->lColElement(i))->NumDofPerElementAuxiliary();
  }

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0,eledofs,0,false));
  acoudishdg->AddDofSet(dofsetaux);

  // call fill complete on acoustical discretization
  acoudishdg->FillComplete();

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
  // if ( DRT::Problem::Instance()->SolverParams(linsolvernumber_acou).get<int>("AZREUSE") < 10 )
  //   dserror("note: you can set AZREUSE to NUMSTEP (ACOUSTIC DYNAMIC) because the system matrix does not change -> save lots of CPU time");

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
    params->set<bool>("acouopt",false);
    params->set<bool>("reduction",false);
  }

  // set restart step if we do not perform inverse analysis
  int restart = DRT::Problem::Instance()->Restart();
  params->set<int>("restart",restart);


  // the main part of this function:
  if(!invanalysis)
  {
    // do we do a PHOTOacoustic simulation? if so -> do a SCATRA simulation for initial pressure distribution
    bool photoacou = DRT::INPUT::IntegralValue<bool>(acouparams,"PHOTOACOU");

    // create algorithm depending on the temporal discretization
    Teuchos::RCP<ACOU::AcouTimeInt> acoualgo;

    INPAR::ACOU::DynamicType dyna = DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(acouparams,"TIMEINT");

    switch(dyna)
    {
    case INPAR::ACOU::acou_impleuler:
    {
      acoualgo = Teuchos::rcp(new ACOU::TimIntImplEuler(acoudishdg,solver,params,output));
      break;
    }
    case INPAR::ACOU::acou_expleuler:
    case INPAR::ACOU::acou_classrk4:
    case INPAR::ACOU::acou_lsrk45reg2:
    case INPAR::ACOU::acou_lsrk33reg2:
    case INPAR::ACOU::acou_lsrk45reg3:
    case INPAR::ACOU::acou_ssprk:
    case INPAR::ACOU::acou_ader:
    case INPAR::ACOU::acou_ader_lts:
    {
      acoualgo = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acoudishdg,solver,params,output));
      break;
    }
    case INPAR::ACOU::acou_ader_tritet:
    {
      acoualgo = Teuchos::rcp(new ACOU::AcouTimeIntAderTriTet(acoudishdg,solver,params,output));
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
      if(0)
      {
        redistributescatradis(scatradis,acoudishdg);
        scatradis->FillComplete();
      }

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

          // add proxy of velocity related degrees of freedom to scatra discretization
          Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
              Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim()+1, 0, 0, true));
          if ( scatradis->AddDofSet(dofsetaux)!= 1 )
            dserror("Scatra discretization has illegal number of dofsets!");

          // finalize discretization
          scatradis->FillComplete(true, false, false);

          // get linear solver id from SCALAR TRANSPORT DYNAMIC
          const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
          if (linsolvernumber == (-1))
            dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

          // create instance of scalar transport basis algorithm (empty fluid discretization)
          Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly =
              Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

          // now we can call Init() on the base algorithm
          // time integrator is constructed and initialized inside
          scatraonly->Init(
              scatradyn,
              scatradyn,
              DRT::Problem::Instance()->SolverParams(linsolvernumber));

          // now me may redistribute or ghost the scatra discretization

          // only now we must call Setup() on th base algorithm
          // all objects relying on the parallel distribution are
          // created and pointers are set.
          // calls Setup() on time integrator inside.
          scatraonly->Setup();

          // set velocity field
          //(this is done only once. Time-dependent velocity fields are not supported)
          (scatraonly->ScaTraField())->SetVelocityField(1);

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

    acoualgo->PrintInformationToScreen();

    // set initial field
    if (restart) // standard restart scenario
    {
      acoualgo->ReadRestart(restart);
    }
    else if (!photoacou) // standard acoustic simulation
    {
      int startfuncno = acouparams.get<int>("STARTFUNCNO");
      acoualgo->SetInitialField(startfuncno);
    }
    else // initial field calculated in dependence on light distribution - photoacoustic setting
    {
      bool meshconform = DRT::INPUT::IntegralValue<bool>(acouparams,"MESHCONFORM");
      // calculate initial pressure distribution by means of the scalar transport problem
      acoualgo->SetInitialPhotoAcousticField(scatrainitialpress,DRT::Problem::Instance()->GetDis("scatra"), meshconform);
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
    Teuchos::RCP<Teuchos::ParameterList> scatraparams = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ScalarTransportDynamicParams()));

    // access the scatra discretization
    Teuchos::RCP<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");

    // ensure that all dofs are assigned in the right order
    scatradis->FillComplete();

    //if ( scatradis->NumGlobalElements() == 0 )
    //  dserror("you said you want to do photoacoustics but you did not supply TRANSP elements");

    // add proxy of velocity related degrees of freedom to scatra discretization
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim()+1, 0, 0, true));
    if ( scatradis->AddDofSet(dofsetaux)!= 1 )
      dserror("Scatra discretization has illegal number of dofsets!");

    // finalize discretization
    scatradis->FillComplete(true, false, false);

    // till here it was all the same...
    // now, we set up the inverse analysis algorithm
    const int linsolvernumber_scatra = scatraparams->get<int>("LINEAR_SOLVER");
    if ( linsolvernumber_scatra == (-1) )
      dserror("no linear solver defined for acoustical problem. Please set LINEAR_SOLVER in ACOUSTIC DYNAMIC to a valid number!");

    // create solver
    Teuchos::RCP<LINALG::Solver> scatrasolver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber_scatra),
                                      scatradis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));

    // create scatra output
    Teuchos::RCP<IO::DiscretizationWriter> scatraoutput = scatradis->Writer();

    // create and run the inverse problem
    Teuchos::RCP<ACOU::PatImageReconstruction> myinverseproblem;
    switch(DRT::INPUT::IntegralValue<INPAR::ACOU::InvAnalysisType>(acouparams,"INV_ANALYSIS"))
    {
    case INPAR::ACOU::pat_none:
      dserror("you should not be here");
    break;
    case INPAR::ACOU::pat_opti:
      myinverseproblem = Teuchos::rcp(new ACOU::PatImageReconstruction(scatradis,acoudishdg,scatraparams,params,scatrasolver,solver,scatraoutput,output));
    break;
    case INPAR::ACOU::pat_optisplit:
      myinverseproblem = Teuchos::rcp(new ACOU::PatImageReconstructionOptiSplit(scatradis,acoudishdg,scatraparams,params,scatrasolver,solver,scatraoutput,output));
    break;
    case INPAR::ACOU::pat_optisplitacousplit:
      myinverseproblem = Teuchos::rcp(new ACOU::PatImageReconstructionOptiSplitAcouSplit(scatradis,acoudishdg,scatraparams,params,scatrasolver,solver,scatraoutput,output));
    break;
    case INPAR::ACOU::pat_optisplitacouident:
      myinverseproblem = Teuchos::rcp(new ACOU::PatImageReconstructionOptiSplitAcouIdent(scatradis,acoudishdg,scatraparams,params,scatrasolver,solver,scatraoutput,output));
    break;
    case INPAR::ACOU::pat_reduction:
      myinverseproblem = Teuchos::rcp(new ACOU::PatImageReconstructionReduction(scatradis,acoudishdg,scatraparams,params,scatrasolver,solver,scatraoutput,output));
    break;
    default:
      dserror("other pat types not listed");
      break;
    }
    if(DRT::INPUT::IntegralValue<bool>(acouparams.sublist("PA IMAGE RECONSTRUCTION"),"SAMPLEOBJECTIVE"))
      myinverseproblem->SampleObjectiveFunction();

    if(restart>0)
      myinverseproblem->ReadRestart(restart);

    myinverseproblem->Optimize();

    // do testing if required
    DRT::Problem::Instance()->AddFieldTest(myinverseproblem->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(scatradis->Comm());

  } // else ** if(!invanalysis)
  return;

}

