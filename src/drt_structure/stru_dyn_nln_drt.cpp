/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "strugenalpha.H"
#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_structure.H"
#include "stru_resulttest.H"

#include "str_invanalysis.H"
#include "../drt_inv_analysis/inv_analysis.H"

#include "../drt_statmech/statmech_time.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
extern "C"
void caldyn_drt()
{
  // get input lists
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  // do we want to do inverse analysis?
  if ((bool)Teuchos::getIntegralValue<int>(iap,"INV_ANALYSIS"))
  {
    STR::invanalysis();
  }
  else
  {
    // get input lists
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    // major switch to different time integrators
    switch (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
    {
    case INPAR::STR::dyna_centr_diff:
      dserror("no central differences in DRT");
      break;
    case INPAR::STR::dyna_gen_alfa:
    case INPAR::STR::dyna_gen_alfa_statics:
      dyn_nlnstructural_drt();
      break;
    case INPAR::STR::dyna_Gen_EMM:
      dserror("GEMM not supported");
      break;
    case INPAR::STR::dyna_statics:
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma :
    case INPAR::STR::dyna_euimsto :
      // direct time integration
      STR::strudyn_direct();
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
    }
  }
}


/*----------------------------------------------------------------------*
  | structural nonlinear dynamics (gen-alpha)              m.gee 12/06  |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(actdis);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->StructuralContactParams();
  const Teuchos::ParameterList& statmech = DRT::Problem::Instance()->StatisticalMechanicsParams();
  const Teuchos::ParameterList& iap      = DRT::Problem::Instance()->InverseAnalysisParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  LINALG::Solver solver(DRT::Problem::Instance()->StructSolverParams(),
                        actdis->Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle());
  actdis->ComputeNullSpaceIfNecessary(solver.Params());

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  switch (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
    //==================================================================
    // Generalized alpha time integration
    //==================================================================
    case INPAR::STR::dyna_gen_alfa :
    case INPAR::STR::dyna_gen_alfa_statics :
    {
      ParameterList genalphaparams;
      StruGenAlpha::SetDefaults(genalphaparams);

      genalphaparams.set<string>("DYNAMICTYP",sdyn.get<string>("DYNAMICTYP"));

      INPAR::STR::ControlType controltype = Teuchos::getIntegralValue<INPAR::STR::ControlType>(sdyn,"CONTROLTYPE");
      genalphaparams.set<INPAR::STR::ControlType>("CONTROLTYPE",controltype);

      {
        vector<int> controlnode;
        std::istringstream contnode(Teuchos::getNumericStringParameter(sdyn,"CONTROLNODE"));
        std::string word;
        while (contnode >> word)
          controlnode.push_back(std::atoi(word.c_str()));
        if ((int)controlnode.size() != 3) dserror("Give proper values for CONTROLNODE in input file");
        genalphaparams.set("CONTROLNODE",controlnode[0]);
        genalphaparams.set("CONTROLDOF",controlnode[1]);
        genalphaparams.set("CONTROLCURVE",controlnode[2]);
      }

      {
        // use linearization of follower loads in Newton
        int loadlin = Teuchos::getIntegralValue<int>(sdyn,"LOADLIN");
        genalphaparams.set<bool>("LOADLIN",loadlin!=0);
      }

      // Rayleigh damping
      genalphaparams.set<bool>  ("damping",(not (sdyn.get<std::string>("DAMPING") == "no"
                                                 or sdyn.get<std::string>("DAMPING") == "No"
                                                 or sdyn.get<std::string>("DAMPING") == "NO")));
      genalphaparams.set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
      genalphaparams.set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

      // Generalised-alpha coefficients
      genalphaparams.set<double>("beta",sdyn.get<double>("BETA"));
#ifdef STRUGENALPHA_BE
      genalphaparams.set<double>("delta",sdyn.get<double>("DELTA"));
#endif
      genalphaparams.set<double>("gamma",sdyn.get<double>("GAMMA"));
      genalphaparams.set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
      genalphaparams.set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

      genalphaparams.set<double>("total time",0.0);
      genalphaparams.set<double>("delta time",sdyn.get<double>("TIMESTEP"));
      genalphaparams.set<double>("max time",sdyn.get<double>("MAXTIME"));
      genalphaparams.set<int>   ("step",0);
      genalphaparams.set<int>   ("nstep",sdyn.get<int>("NUMSTEP"));
      genalphaparams.set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
      genalphaparams.set<int>   ("num iterations",-1);

      genalphaparams.set<string>("convcheck", sdyn.get<string>("CONV_CHECK"));
      genalphaparams.set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
      genalphaparams.set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
      genalphaparams.set<double>("tolerance constraint",sdyn.get<double>("TOLCONSTR"));

      genalphaparams.set<double>("UZAWAPARAM",sdyn.get<double>("UZAWAPARAM"));
      genalphaparams.set<double>("UZAWATOL",sdyn.get<double>("UZAWATOL"));
      genalphaparams.set<int>   ("UZAWAMAXITER",sdyn.get<int>("UZAWAMAXITER"));
      genalphaparams.set<INPAR::STR::ConSolveAlgo>("UZAWAALGO",getIntegralValue<INPAR::STR::ConSolveAlgo>(sdyn,"UZAWAALGO"));
      genalphaparams.set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
      genalphaparams.set<int>   ("io disp every nstep",sdyn.get<int>("RESEVRYDISP"));

      genalphaparams.set<bool>  ("ADAPTCONV",getIntegralValue<int>(sdyn,"ADAPTCONV")==1);
      genalphaparams.set<double>("ADAPTCONV_BETTER",sdyn.get<double>("ADAPTCONV_BETTER"));

      INPAR::STR::StressType iostress = Teuchos::getIntegralValue<INPAR::STR::StressType>(ioflags,"STRUCT_STRESS");
      genalphaparams.set<INPAR::STR::StressType>("io structural stress", iostress);

      genalphaparams.set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

      INPAR::STR::StrainType iostrain = Teuchos::getIntegralValue<INPAR::STR::StrainType>(ioflags,"STRUCT_STRAIN");
      genalphaparams.set<INPAR::STR::StrainType>("io structural strain", iostrain);

      genalphaparams.set<bool>  ("io surfactant",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_SURFACTANT"));

      genalphaparams.set<int>   ("restart",probtype.get<int>("RESTART"));
      genalphaparams.set<int>   ("write restart every",sdyn.get<int>("RESTARTEVRY"));

      genalphaparams.set<bool>  ("print to screen",true);
      genalphaparams.set<bool>  ("print to err",true);
      genalphaparams.set<FILE*> ("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

      // parameters for inverse analysis
      genalphaparams.set<bool>  ("inv_analysis",Teuchos::getIntegralValue<int>(iap,"INV_ANALYSIS"));
      genalphaparams.set<double>("measured_curve0",iap.get<double>("MEASURED_CURVE0"));
      genalphaparams.set<double>("measured_curve1",iap.get<double>("MEASURED_CURVE1"));
      genalphaparams.set<double>("measured_curve2",iap.get<double>("MEASURED_CURVE2"));
      genalphaparams.set<double>("inv_ana_tol",iap.get<double>("INV_ANA_TOL"));
      
      // non-linear solution technique
      switch (Teuchos::getIntegralValue<INPAR::STR::NonlinSolTech>(sdyn,"NLNSOL"))
      {
        case INPAR::STR::soltech_newtonfull:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
        case INPAR::STR::soltech_newtonls:
          genalphaparams.set<string>("equilibrium iteration","line search newton");
        break;
        case INPAR::STR::soltech_newtonopp:
          genalphaparams.set<string>("equilibrium iteration","oppositely converging newton");
        break;
        case INPAR::STR::soltech_newtonmod:
          genalphaparams.set<string>("equilibrium iteration","modified newton");
        break;
        case INPAR::STR::soltech_nlncg:
          genalphaparams.set<string>("equilibrium iteration","nonlinear cg");
        break;
        case INPAR::STR::soltech_ptc:
          genalphaparams.set<string>("equilibrium iteration","ptc");
        break;
        case INPAR::STR::soltech_newtonuzawalin:
          genalphaparams.set<string>("equilibrium iteration","newtonlinuzawa");
        break;
        case INPAR::STR::soltech_newtonuzawanonlin:
          genalphaparams.set<string>("equilibrium iteration","augmentedlagrange");
        break;
        default:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
      }

      // set predictor (takes values "constant" "consistent")
      switch (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT"))
      {
        case INPAR::STR::pred_vague:
          dserror("You have to define the predictor");
          break;
        case INPAR::STR::pred_constdis:
          genalphaparams.set<string>("predictor","consistent");
          break;
        case INPAR::STR::pred_constdisvelacc:
          genalphaparams.set<string>("predictor","constant");
          break;
        case INPAR::STR::pred_tangdis:
          genalphaparams.set<string>("predictor","tangdis");
          break;
        default:
          dserror("Cannot cope with choice of predictor");
          break;
      }

      // detect if contact is present
      bool contact = false;
      switch (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(scontact,"CONTACT"))
      {
        case INPAR::CONTACT::contact_none:
          contact = false;
          break;
        case INPAR::CONTACT::contact_normal:
          contact = true;
          break;
        case INPAR::CONTACT::contact_frictional:
          contact = true;
          break;
        case INPAR::CONTACT::contact_meshtying:
          contact = true;
          break;
        default:
          dserror("Cannot cope with choice of contact type");
          break;
      }

      // detect whether thermal bath is present
      bool thermalbath = false;
      switch (Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmech,"THERMALBATH"))
      {
        case INPAR::STATMECH::thermalbath_none:
          thermalbath = false;
          break;
        case INPAR::STATMECH::thermalbath_uniform:
          thermalbath = true;
          break;
        case INPAR::STATMECH::thermalbath_shearflow:
          thermalbath = true;
          break;
        default:
          dserror("Cannot cope with choice of thermal bath");
          break;
      }

      // create the time integrator
      bool inv_analysis = genalphaparams.get("inv_analysis",false);
      RCP<StruGenAlpha> tintegrator;
      if (!contact && !inv_analysis && !thermalbath)
        tintegrator = rcp(new StruGenAlpha(genalphaparams,*actdis,solver,output));
      else
      {
        if (contact)
          tintegrator = rcp(new CONTACT::ContactStruGenAlpha(genalphaparams,*actdis,solver,output));
        if (inv_analysis)
          dserror("Inverse analysis moved ahead to STI");
        if (thermalbath)
          tintegrator = rcp(new StatMechTime(genalphaparams,*actdis,solver,output));
      }

      // do restart if demanded from input file
      // note that this changes time and step in genalphaparams
      if (genprob.restart)
        tintegrator->ReadRestart(genprob.restart);

      // write mesh always at beginning of calc or restart
      {
        int    step = genalphaparams.get<int>("step",0);
        double time = genalphaparams.get<double>("total time",0.0);
        output.WriteMesh(step,time);
      }

      // integrate in time and space
      tintegrator->Integrate();

      // test results
      {
        DRT::Problem::Instance()->AddFieldTest(rcp(new StruResultTest(*tintegrator)));
        DRT::Problem::Instance()->TestAll(actdis->Comm());
      }

    }
    break;
    //==================================================================
    // Generalized Energy Momentum Method
    //==================================================================
    case INPAR::STR::dyna_Gen_EMM :
    {
      dserror("Not yet impl.");
    }
    break;
    //==================================================================
    // Everything else
    //==================================================================
    default :
    {
      dserror("Time integration scheme is not available");
    }
    break;
  } // end of switch(sdyn->Typ)

  Teuchos::TimeMonitor::summarize();

  return;
} // end of dyn_nlnstructural_drt()

#endif  // #ifdef CCADISCRET
