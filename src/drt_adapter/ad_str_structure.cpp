/*----------------------------------------------------------------------*/
/*!
\file adapter_structure.cpp

\brief Structure field adapter

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_str_structure.H"
#include "ad_str_timint_adaptive.H"
#include "ad_str_constr_merged.H"
#include "ad_str_wrapper.H"
#include "ad_str_lung.H"
#include "ad_str_redairway.H"
#include "ad_str_fsi_crack.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_meshfree.H"
#include "../drt_inpar/inpar_crack.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_io/io_control.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimada_create.H"
#include "../drt_structure/strtimint_impl.H"
#include "../drt_patspec/patspec.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"

#include "../drt_io/io_pstream.H"

#include "../drt_crack/InsertCohesiveElements.H"
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Structure::~Structure()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::StructureBaseAlgorithm(
  const Teuchos::ParameterList& prbdyn,
  Teuchos::ParameterList& sdyn,
  Teuchos::RCP<DRT::Discretization> actdis
)
{
  SetupStructure(prbdyn, sdyn, actdis);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::~StructureBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupStructure(
  const Teuchos::ParameterList& prbdyn,
  Teuchos::ParameterList& sdyn,
  Teuchos::RCP<DRT::Discretization> actdis
)
{
  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
  case INPAR::STR::dyna_statics :
  case INPAR::STR::dyna_genalpha :
  case INPAR::STR::dyna_onesteptheta :
  case INPAR::STR::dyna_gemm :
  case INPAR::STR::dyna_expleuler:
  case INPAR::STR::dyna_centrdiff :
  case INPAR::STR::dyna_particle_centrdiff :
  case INPAR::STR::dyna_ab2 :
  case INPAR::STR::dyna_euma :
  case INPAR::STR::dyna_euimsto :
  case INPAR::STR::dyna_statmech :
    SetupTimInt(prbdyn, sdyn, actdis);  // <-- here is the show
    break;
  default :
    dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
    break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupTimInt(
  const Teuchos::ParameterList& prbdyn,
  Teuchos::ParameterList& sdyn,
  Teuchos::RCP<DRT::Discretization> actdis
)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t
    = Teuchos::TimeMonitor::getNewTimer("ADAPTER::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // get the problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  // what's the current problem type?
  PROBLEM_TYP probtype = problem->ProblemType();

  // set degrees of freedom in the discretization
  if (not actdis->Filled() || not actdis->HaveDofs()) actdis->FillComplete();

  // get input parameter lists and copy them, because a few parameters are overwritten
  //const Teuchos::ParameterList& probtype
  //  = problem->ProblemTypeParams();
  Teuchos::RCP<Teuchos::ParameterList> ioflags
    = Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<Teuchos::ParameterList> tap
    = Teuchos::rcp(new Teuchos::ParameterList(sdyn.sublist("TIMEADAPTIVITY")));
  Teuchos::RCP<Teuchos::ParameterList> snox
    = Teuchos::rcp(new Teuchos::ParameterList(problem->StructuralNoxParams()));

  // show default parameters
  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", problem->ErrorFile()->Handle());
  Teuchos::ParameterList& nox = xparams->sublist("NOX");
  nox = *snox;
  // Parameter to determine if MLMC is on/off
  Teuchos::RCP<Teuchos::ParameterList> mlmcp
      = Teuchos::rcp(new Teuchos::ParameterList (problem->MultiLevelMonteCarloParams()));
  // Needed for reduced restart output
  xparams->set<int>("REDUCED_OUTPUT",Teuchos::getIntegralValue<int>((*mlmcp),"REDUCED_OUTPUT"));

  //dual problems need the primal results which need to be stored therefore
  const Teuchos::ParameterList& ivap = problem->StatInverseAnalysisParams();
  xparams->set<int>("MSTEPEVRY",Teuchos::getIntegralValue<int>(ivap,"MSTEPS"));


  // overrule certain parameters
  sdyn.set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  sdyn.set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  sdyn.set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  if(probtype == prb_struct_ale || probtype == prb_structure || probtype == prb_redairways_tissue || probtype == prb_particle)
  {
    sdyn.set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));
  }
  else
  {
    sdyn.set<int>("RESULTSEVRY", prbdyn.get<int>("UPRES"));
  }

  // Check if for chosen Rayleigh damping the regarding parameters are given explicitly in the .dat file
  if (DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(sdyn,"DAMPING") == INPAR::STR::damp_rayleigh)
    {
  		if (sdyn.get<double>("K_DAMP") < 0.0)
  		  {
  		    dserror("Rayleigh damping parameter K_DAMP not explicitly given.");
  		  }
  		if (sdyn.get<double>("M_DAMP") < 0.0)
      {
         dserror("Rayleigh damping parameter M_DAMP not explicitly given.");
       }
    }

  // create a solver
  Teuchos::RCP<LINALG::Solver> solver = CreateLinearSolver(actdis, sdyn);

  // create contact/meshtying solver only if contact/meshtying problem.
  Teuchos::RCP<LINALG::Solver> contactsolver = Teuchos::null;
  const Teuchos::ParameterList& scontact = problem->ContactDynamicParams();
  INPAR::CONTACT::ApplicationType bContact = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION");
  switch (bContact)
  {
    case INPAR::CONTACT::app_none:
      break;
    case INPAR::CONTACT::app_mortarcontact:
    case INPAR::CONTACT::app_mortarmeshtying:
    case INPAR::CONTACT::app_beamcontact:
      contactsolver = CreateContactMeshtyingSolver(actdis, sdyn);
      break;
    default:
      dserror("Cannot cope with choice of contact or meshtying type");
      break;
  }

  // create contact/meshtying solver only if contact/meshtying problem.
  //Teuchos::RCP<LINALG::Solver> contactsolver = CreateContactMeshtyingSolver(actdis, sdyn);

  if (solver != Teuchos::null &&
      (solver->Params().isSublist("Aztec Parameters") || solver->Params().isSublist("Belos Parameters"))
      &&
      solver->Params().isSublist("ML Parameters") // TODO what about MueLu?
      &&
      DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdyn,"STC_SCALING")!=INPAR::STR::stc_none)
    {
      Teuchos::ParameterList& mllist = solver->Params().sublist("ML Parameters");
      Teuchos::RCP<std::vector<double> > ns = mllist.get<Teuchos::RCP<std::vector<double> > >("nullspace");

      const int size=actdis->DofRowMap()->NumMyElements();

      // extract the six nullspace vectors corresponding to the modes
      // trans x, trans y, trans z, rot x, rot y, rot z
      // Note: We assume 3d here!

      Teuchos::RCP<Epetra_Vector> nsv1=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[0])));
      Teuchos::RCP<Epetra_Vector> nsv2=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[size])));
      Teuchos::RCP<Epetra_Vector> nsv3=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[2*size])));
      Teuchos::RCP<Epetra_Vector> nsv4=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[3*size])));
      Teuchos::RCP<Epetra_Vector> nsv5=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[4*size])));
      Teuchos::RCP<Epetra_Vector> nsv6=
          Teuchos::rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[5*size])));


      //prepare matrix for scaled thickness business of thin shell structures
      Teuchos::RCP<LINALG::SparseMatrix> stcinv=
        Teuchos::rcp(new LINALG::SparseMatrix(*actdis->DofRowMap(), 81, true, true));

      stcinv->Zero();
      // create the parameters for the discretization
      Teuchos::ParameterList p;
      // action for elements
      const std::string action = "calc_stc_matrix_inverse";
      p.set("action", action);
      p.set<int>("stc_scaling",DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdyn,"STC_SCALING"));
      p.set("stc_layer",1);

      actdis-> Evaluate(p, stcinv, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);

      stcinv->Complete();

      for (int lay = 2; lay <= sdyn.get<int>("STC_LAYER"); ++lay)
      {
        Teuchos::ParameterList pe;

        p.set("stc_layer", lay);

        Teuchos::RCP<LINALG::SparseMatrix> tmpstcmat=
          Teuchos::rcp(new LINALG::SparseMatrix(*actdis->DofRowMap(),81,true,true));
        tmpstcmat->Zero();

        actdis-> Evaluate(p, tmpstcmat, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);
        tmpstcmat->Complete();

        stcinv = MLMultiply(*stcinv,*tmpstcmat,false,false,true);
      }

      Teuchos::RCP<Epetra_Vector> temp =
          LINALG::CreateVector(*(actdis->DofRowMap()),false);

      stcinv->Multiply(false,*nsv1,*temp);
      nsv1->Update(1.0,*temp,0.0);
      stcinv->Multiply(false,*nsv2,*temp);
      nsv2->Update(1.0,*temp,0.0);
      stcinv->Multiply(false,*nsv3,*temp);
      nsv3->Update(1.0,*temp,0.0);
      stcinv->Multiply(false,*nsv4,*temp);
      nsv4->Update(1.0,*temp,0.0);
      stcinv->Multiply(false,*nsv5,*temp);
      nsv5->Update(1.0,*temp,0.0);
      stcinv->Multiply(false,*nsv6,*temp);
      nsv6->Update(1.0,*temp,0.0);

    }

  // Checks in case of multi-scale simulations

  {
    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<MAT::PAR::Bundle> materials = problem->Materials();
    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials->Map()->begin();
         i!=materials->Map()->end();
         ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> mat = i->second;
      if (mat->Type() == INPAR::MAT::m_struct_multiscale)
      {
        if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP") != INPAR::STR::dyna_genalpha)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(sdyn.sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_trlike)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=TrLike");
        break;
      }
    }
  }

  // Add cohesive elements in case of crack propagation simulations
  {
    const Teuchos::ParameterList& crackparam = DRT::Problem::Instance()->CrackParams();
    if (DRT::INPUT::IntegralValue<INPAR::CRACK::CrackModel>(crackparam,"CRACK_MODEL")
          != INPAR::CRACK::crack_none)
    {
      //if( (not DRT::Problem::Instance()->Restart()) and DRT::Problem::Instance()->ProblemType() == prb_structure )
      {
        Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
        DRT::CRACK::InsertCohesiveElements isp( structdis );
      }
    }
  }

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  bool writeinitialmesh = DRT::INPUT::IntegralValue<bool>(sdyn,"WRITE_INITIAL_MESH");
  if (writeinitialmesh)
  {
    output->WriteMesh(0, 0.0);
  }

  // create marching time integrator
  Teuchos::RCP<STR::TimInt> tmpstr = STR::TimIntCreate(*ioflags, sdyn, *xparams, actdis, solver, contactsolver, output);

  // create auxiliar time integrator, can be seen as a wrapper for tmpstr
  Teuchos::RCP<STR::TimAda> sta = STR::TimAdaCreate(*ioflags, sdyn, *xparams, *tap, tmpstr);

  if (sta!=Teuchos::null)
  {
    if (probtype != prb_structure)
    {
      dserror("adaptive time integration for the structure only tested for problem type prb_structure");
    }
    structure_ = Teuchos::rcp(new StructureTimIntAda(sta, tmpstr));
  }
  else if (tmpstr != Teuchos::null)
  {
    switch(probtype)
    {
    case prb_fsi:
    case prb_fsi_lung:
    case prb_gas_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_fsi_xfem:
    case prb_fluid_fluid_fsi:
    {
      const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

      if ((actdis->Comm()).MyPID()==0) IO::cout << "Using StructureNOXCorrectionWrapper()..." << IO::endl;

      if (tmpstr->HaveConstraint())
      {
        if (coupling == fsi_iter_constr_monolithicstructuresplit or
            coupling == fsi_iter_constr_monolithicfluidsplit)
          structure_ = Teuchos::rcp(new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        else
          structure_ = Teuchos::rcp(new StructureConstrMerged(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
      }
      else
      {
        if (coupling == fsi_iter_lung_monolithicstructuresplit or
            coupling == fsi_iter_lung_monolithicfluidsplit)
          structure_ = Teuchos::rcp(new StructureLung(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        else
          structure_ = Teuchos::rcp(new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
      }
    }
    break;
    case prb_fsi_crack:
    {
      structure_ = Teuchos::rcp(new FSICrackingStructure(Teuchos::rcp(new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))))));
    }
    break;
    case prb_redairways_tissue:
    {
      structure_ = Teuchos::rcp(new StructureRedAirway(tmpstr));
    }
    break;
    case prb_poroelast:
    case prb_poroscatra:
    case prb_fpsi:
    {
      const Teuchos::ParameterList& porodyn = problem->PoroelastDynamicParams();
      const INPAR::POROELAST::SolutionSchemeOverFields coupling =
            DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(porodyn, "COUPALGO");
      if (tmpstr->HaveConstraint())
      {
        if (coupling == INPAR::POROELAST::Monolithic_structuresplit or coupling == INPAR::POROELAST::Monolithic_fluidsplit)
          structure_ = Teuchos::rcp(new FSIStructureWrapper(tmpstr));
        else
          structure_ = Teuchos::rcp(new StructureConstrMerged(tmpstr));
      }
      else
      {
          structure_ = Teuchos::rcp(new FSIStructureWrapper(tmpstr));
      }
    }
    break;
    case prb_struct_ale:
    {
      structure_ = Teuchos::rcp(new FSIStructureWrapper(tmpstr));
    }
    break;
    default:
    {
      structure_ = tmpstr;
    }
    break;
    }
  }
  else
  {
    dserror("no proper time integration found");
  }
  // see you
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> ADAPTER::StructureBaseAlgorithm::CreateLinearSolver(Teuchos::RCP<DRT::Discretization>& actdis, const Teuchos::ParameterList& sdyn)
{
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  // get the solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");
  if (linsolvernumber == (-2))
    return Teuchos::null;

  solver = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                    actdis->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  return solver;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> ADAPTER::StructureBaseAlgorithm::CreateContactMeshtyingSolver(Teuchos::RCP<DRT::Discretization>& actdis, const Teuchos::ParameterList& sdyn)
{
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  const Teuchos::ParameterList& mcparams     = DRT::Problem::Instance()->ContactDynamicParams();
//  INPAR::CONTACT::ApplicationType apptype    = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(mcparams,"APPLICATION");
  switch(DRT::INPUT::IntegralValue<int>(mcparams,"SYSTEM"))
  {
    case INPAR::CONTACT::system_spsimpler:
    {
      // meshtying/contact for structure
      // get the solver number used for meshtying/contact problems
      const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
      // check if the meshtying/contact solver has a valid solver number
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // plausibility check
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
      if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
          prec != INPAR::SOLVER::azprec_TekoSIMPLE  &&
          prec != INPAR::SOLVER::azprec_MueLuAMG_contactSP)  // TODO adapt error message
        dserror("Mortar/Contact with saddlepoint system only possible with SIMPLE preconditioner. Choose CheapSIMPLE or TekoSIMPLE in the SOLVER %i block in your dat file.",linsolvernumber);

      // build meshtying solver
      solver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                             actdis->Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

      actdis->ComputeNullSpaceIfNecessary(solver->Params());

      INPAR::CONTACT::ApplicationType apptype = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(mcparams,"APPLICATION");
      if     (apptype == INPAR::CONTACT::app_mortarcontact) solver->Params().set<bool>("CONTACT",true);
      else if(apptype==INPAR::CONTACT::app_mortarmeshtying) solver->Params().set<bool>("MESHTYING",true);
      else dserror("this cannot be: no saddlepoint problem for beamcontact or pure structure problem.");

      INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(mcparams,"STRATEGY");
      INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(mcparams,"SYSTEM");
      if (soltype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed)
      {
        // get the solver number used for structural problems
        const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
        // check if the structural solver has a valid solver number
        if (linsolvernumber == (-1))
          dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

        // provide null space information
        if (prec == INPAR::SOLVER::azprec_CheapSIMPLE ||
                  prec == INPAR::SOLVER::azprec_TekoSIMPLE) {
          actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1")); // Inverse2 is created within blockpreconditioners.cpp
        } else if (prec == INPAR::SOLVER::azprec_MueLuAMG_contactSP) {
          // note: the null space is definitely too long and wrong for the Lagrange multipliers
          // don't forget to call FixMLNullspace for "Inverse1"!
          //actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1")); // TODO do we need this?
        }


      }
    }
    break;
    default:
    {
      // meshtying/contact for structure
       // get the solver number used for meshtying/contact problems
       const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
       // check if the meshtying/contact solver has a valid solver number
       if (linsolvernumber == (-1))
         dserror("no linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // build meshtying solver
      solver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                             actdis->Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));
      actdis->ComputeNullSpaceIfNecessary(solver->Params());
    }
    break;
  }

  return solver;
}

/*----------------------------------------------------------------------*/
