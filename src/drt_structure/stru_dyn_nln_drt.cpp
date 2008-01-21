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

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "stru_genalpha_zienxie_drt.H"
#include "strugenalpha.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../io/io_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
extern "C"
void caldyn_drt()
{
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();

  switch (Teuchos::getIntegralValue<int>(sdyn,"DYNAMICTYP"))
  {
  case STRUCT_DYNAMIC::centr_diff:
    dserror("no central differences in DRT");
    break;
  case STRUCT_DYNAMIC::gen_alfa:
    switch (Teuchos::getIntegralValue<int>(sdyn,"TA_KIND"))
    {
    case TIMADA_DYNAMIC::timada_kind_none:
      dyn_nlnstructural_drt();
      break;
    case TIMADA_DYNAMIC::timada_kind_zienxie:
      stru_genalpha_zienxie_drt();
      break;
    default:
      dserror("unknown time adaption scheme '%s'", sdyn.get<std::string>("TA_KIND").c_str());
    }
    break;
  case STRUCT_DYNAMIC::Gen_EMM:
    dserror("GEMM not supported");
    break;
  default:
    dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
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
  SOLVAR*         actsolv  = &solv[0];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  switch (Teuchos::getIntegralValue<int>(sdyn,"DYNAMICTYP"))
  {
    //==================================================================
    // Generalized alpha time integration
    //==================================================================
    case STRUCT_DYNAMIC::gen_alfa :
    {
      ParameterList genalphaparams;
      StruGenAlpha::SetDefaults(genalphaparams);

      genalphaparams.set<bool>  ("damping",Teuchos::getIntegralValue<int>(sdyn,"DAMPING"));
      genalphaparams.set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
      genalphaparams.set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

      genalphaparams.set<double>("beta",sdyn.get<double>("BETA"));
      genalphaparams.set<double>("gamma",sdyn.get<double>("GAMMA"));
      genalphaparams.set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
      genalphaparams.set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

      genalphaparams.set<double>("total time",0.0);
      genalphaparams.set<double>("delta time",sdyn.get<double>("TIMESTEP"));
      genalphaparams.set<int>   ("step",0);
      genalphaparams.set<int>   ("nstep",sdyn.get<int>("NUMSTEP"));
      genalphaparams.set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
      genalphaparams.set<int>   ("num iterations",-1);

      genalphaparams.set<string>("convcheck", sdyn.get<string>("CONV_CHECK"));
      genalphaparams.set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
      genalphaparams.set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
      genalphaparams.set<double>("tolerance volume",sdyn.get<double>("TOLVOL"));

      genalphaparams.set<bool>  ("contact",static_cast<bool>(sdyn.get<int>("CONTACT")));

      genalphaparams.set<double>("uzawa parameter",sdyn.get<double>("UZAWAPARAM"));
      genalphaparams.set<int>("uzawa maxiter",sdyn.get<int>("UZAWAMAXITER"));
      genalphaparams.set<string>("uzawa algorithm",sdyn.get<string>("UZAWAALGO"));
      genalphaparams.set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
      genalphaparams.set<int>   ("io disp every nstep",sdyn.get<int>("RESEVRYDISP"));
      genalphaparams.set<bool>  ("io structural stress",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_STRESS"));
      genalphaparams.set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

      genalphaparams.set<int>   ("restart",probtype.get<int>("RESTART"));
      genalphaparams.set<int>   ("write restart every",sdyn.get<int>("RESTARTEVRY"));

      genalphaparams.set<bool>  ("print to screen",true);
      genalphaparams.set<bool>  ("print to err",true);
      genalphaparams.set<FILE*> ("err file",allfiles.out_err);

      switch (Teuchos::getIntegralValue<int>(sdyn,"NLNSOL"))
      {
        case STRUCT_DYNAMIC::fullnewton:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
        case STRUCT_DYNAMIC::modnewton:
          genalphaparams.set<string>("equilibrium iteration","modified newton");
        break;
        case STRUCT_DYNAMIC::matfreenewton:
          genalphaparams.set<string>("equilibrium iteration","matrixfree newton");
        break;
        case STRUCT_DYNAMIC::nlncg:
          genalphaparams.set<string>("equilibrium iteration","nonlinear cg");
        break;
        case STRUCT_DYNAMIC::ptc:
          genalphaparams.set<string>("equilibrium iteration","ptc");
        break;
        default:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
      }

      // set predictor (takes values "constant" "consistent")
      switch (Teuchos::getIntegralValue<int>(sdyn,"PREDICT"))
      {
        case STRUCT_DYNAMIC::pred_vague:
          dserror("You have to define the predictor");
          break;
        case STRUCT_DYNAMIC::pred_constdis:
          genalphaparams.set<string>("predictor","consistent");
          break;
        case STRUCT_DYNAMIC::pred_constdisvelacc:
          genalphaparams.set<string>("predictor","constant");
          break;
        default:
          dserror("Cannot cope with choice of predictor");
          break;
      }

      // create the time integrator
      bool contact = genalphaparams.get("contact",false);
      RCP<StruGenAlpha> tintegrator = null;
      if (!contact) tintegrator = rcp(new StruGenAlpha(genalphaparams,*actdis,solver,output));
      else          tintegrator = rcp(new ContactStruGenAlpha(genalphaparams,*actdis,solver,output));

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
    }
    break;
    //==================================================================
    // Generalized Energy Momentum Method
    //==================================================================
    case STRUCT_DYNAMIC::Gen_EMM :
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

  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef CCADISCRET
