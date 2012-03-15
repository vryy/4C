/*----------------------------------------------------------------------*/
/*!
\file adapter_structure.cpp

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_structure.H"
#include "adapter_structure_timint_adaptive.H"
#include "adapter_structure_constr_merged.H"
#include "adapter_structure_wrapper.H"
#include "adapter_structure_lung.H"
#include "adapter_structure_timint_poroelast.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_io/io_control.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimada_create.H"
#include "../drt_beamcontact/beam3contactstrugenalpha.H"
#include "../drt_statmech/statmech_time.H"
#include "../drt_patspec/patspec.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Structure::~Structure()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::StructureBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  SetupStructure(prbdyn);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::~StructureBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupStructure(const Teuchos::ParameterList& prbdyn)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
  case INPAR::STR::dyna_statics :
  case INPAR::STR::dyna_genalpha :
  case INPAR::STR::dyna_onesteptheta :
  case INPAR::STR::dyna_gemm :
  case INPAR::STR::dyna_expleuler:
  case INPAR::STR::dyna_centrdiff :
  case INPAR::STR::dyna_ab2 :
  case INPAR::STR::dyna_euma :
  case INPAR::STR::dyna_euimsto :
    SetupTimIntImpl(prbdyn);  // <-- here is the show
    break;
  default :
    dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
    break;
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupTimIntImpl(const Teuchos::ParameterList& prbdyn)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t
    = Teuchos::TimeMonitor::getNewTimer("ADAPTER::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (not actdis->Filled() || not actdis->HaveDofs()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  // for multilevel monte carlo we do not need to write mesh in every run
  Teuchos::RCP<Teuchos::ParameterList> mlmcp
      = rcp(new Teuchos::ParameterList (DRT::Problem::Instance()->MultiLevelMonteCarloParams()));
  bool perform_mlmc = Teuchos::getIntegralValue<int>((*mlmcp),"MLMC");
    if (perform_mlmc!=true)
    {
      output->WriteMesh(0, 0.0);
    }
    // get input parameter lists and copy them, because a few parameters are overwritten
  //const Teuchos::ParameterList& probtype
  //  = DRT::Problem::Instance()->ProblemTypeParams();
  Teuchos::RCP<Teuchos::ParameterList> ioflags
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->IOParams()));
  Teuchos::RCP<Teuchos::ParameterList> sdyn
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->StructuralDynamicParams()));
  Teuchos::RCP<Teuchos::ParameterList> tap
    = Teuchos::rcp(new Teuchos::ParameterList(sdyn->sublist("TIMEADAPTIVITY")));
  Teuchos::RCP<Teuchos::ParameterList> snox
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->StructuralNoxParams()));

  // show default parameters
  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, *sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  Teuchos::ParameterList& nox = xparams->sublist("NOX");
  nox = *snox;
  // Parameter to determinen wether MLMC is on/off
  // Needed to for reduced restart output
  xparams->set<int>("REDUCED_OUTPUT",Teuchos::getIntegralValue<int>((*mlmcp),"REDUCED_OUTPUT"));

  // overrule certain parameters
  sdyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  sdyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  sdyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  if(genprob.probtyp == prb_struct_ale || genprob.probtyp == prb_structure)
  {
    sdyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));
  }
  else
  {
    sdyn->set<int>("RESULTSEVRY", prbdyn.get<int>("UPRES"));
  }
  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi or
      genprob.probtyp == prb_fsi_lung or
      genprob.probtyp == prb_gas_fsi or
      genprob.probtyp == prb_biofilm_fsi or
      genprob.probtyp == prb_thermo_fsi or
      genprob.probtyp == prb_fluid_fluid_fsi)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn
      = DRT::Problem::Instance()->FSIDynamicParams();

    // Robin flags
    INPAR::FSI::PartitionedCouplingMethod method
      = DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");
    xparams->set<bool>("structrobin",
                       ( (method==INPAR::FSI::DirichletRobin)
                         or (method==INPAR::FSI::RobinRobin) ));

    // THIS SHOULD GO, OR SHOULDN'T IT?
    xparams->set<double>("alpha s", fsidyn.get<double>("ALPHA_S"));

    // check if predictor fits to FSI algo
    int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if ( (coupling == fsi_iter_monolithicfluidsplit)
         or (coupling == fsi_iter_monolithiclagrange)
         or (coupling == fsi_iter_monolithicstructuresplit)
         or (coupling == fsi_iter_lung_monolithicstructuresplit)
         or (coupling == fsi_iter_lung_monolithicfluidsplit)
         or (coupling == fsi_iter_constr_monolithicfluidsplit)
         or (coupling == fsi_iter_constr_monolithicstructuresplit)
         or (coupling == fsi_iter_mortar_monolithicfluidsplit)
         or (coupling == fsi_iter_mortar_monolithicstructuresplit)
         or (coupling == fsi_iter_fluidfluid_monolithicstructuresplit))
    {
      if ((DRT::INPUT::IntegralValue<INPAR::STR::PredEnum>(*sdyn,"PREDICT")
          != INPAR::STR::pred_constdisvelacc) and
          (DRT::INPUT::IntegralValue<INPAR::STR::PredEnum>(*sdyn,"PREDICT")
          != INPAR::STR::pred_constdisvelaccpres))
      {
        dserror("only constant structure predictor with monolithic FSI possible");
      }
    }
  }

  // extra parameters for poroelasticity
  if (genprob.probtyp == prb_poroelast)
  {
	  //set parameters for poroelasticity
	  sdyn->set<double>("INITPOROSITY", prbdyn.get<double>("INITPOROSITY"));
  }

  // create a solver
  Teuchos::RCP<LINALG::Solver> solver = CreateLinearSolver(actdis);

  // create contact/meshtying solver only if contact/meshtying problem.
  Teuchos::RCP<LINALG::Solver> contactsolver = Teuchos::null;
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->MeshtyingAndContactParams();
  INPAR::CONTACT::ApplicationType bContact = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION");
  switch (bContact)
  {
    case INPAR::CONTACT::app_none:
      break;
    case INPAR::CONTACT::app_mortarcontact:
    case INPAR::CONTACT::app_mortarmeshtying:
    case INPAR::CONTACT::app_beamcontact:
      contactsolver = CreateContactMeshtyingSolver(actdis);
      break;
    default:
      dserror("Cannot cope with choice of contact or meshtying type");
      break;
  }

  // create contact/meshtying solver only if contact/meshtying problem.
  //Teuchos::RCP<LINALG::Solver> contactsolver = CreateContactMeshtyingSolver(actdis);

  if ((solver->Params().isSublist("Aztec Parameters") || solver->Params().isSublist("Belos Parameters"))
      &&
      solver->Params().isSublist("ML Parameters") // TODO what about MueLu?
      &&
      DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(*sdyn,"STC_SCALING")!=INPAR::STR::stc_none)
    {
      ParameterList& mllist = solver->Params().sublist("ML Parameters");
      RCP<vector<double> > ns = mllist.get<RCP<vector<double> > >("nullspace");

      const int size=actdis->DofRowMap()->NumMyElements();

      // extract the six nullspace vectors corresponding to the modes
      // trans x, trans y, trans z, rot x, rot y, rot z
      // Note: We assume 3d here!

      Teuchos::RCP<Epetra_Vector> nsv1=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[0])));
      Teuchos::RCP<Epetra_Vector> nsv2=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[size])));
      Teuchos::RCP<Epetra_Vector> nsv3=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[2*size])));
      Teuchos::RCP<Epetra_Vector> nsv4=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[3*size])));
      Teuchos::RCP<Epetra_Vector> nsv5=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[4*size])));
      Teuchos::RCP<Epetra_Vector> nsv6=
          rcp(new Epetra_Vector(View,*(actdis->DofRowMap()),&((*ns)[5*size])));


      //prepare matrix for scaled thickness business of thin shell structures
      Teuchos::RCP<LINALG::SparseMatrix> stcinv=
        Teuchos::rcp(new LINALG::SparseMatrix(*actdis->DofRowMap(), 81, true, true));

      stcinv->Zero();
      // create the parameters for the discretization
      Teuchos::ParameterList p;
      // action for elements
      const std::string action = "calc_stc_matrix_inverse";
      p.set("action", action);
      p.set<int>("stc_scaling",DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(*sdyn,"STC_SCALING"));
      p.set("stc_layer",1);

      actdis-> Evaluate(p, stcinv, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);

      stcinv->Complete();

      for (int lay = 2; lay <= sdyn->get<int>("STC_LAYER"); ++lay)
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
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<MAT::PAR::Bundle> materials = DRT::Problem::Instance()->Materials();
    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials->Map()->begin();
         i!=materials->Map()->end();
         ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> mat = i->second;
      if (mat->Type() == INPAR::MAT::m_struct_multiscale)
      {
        if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP") != INPAR::STR::dyna_genalpha)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(sdyn.sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_imrlike)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=ImrLike");
        break;
      }
    }
  }

  // create marching time integrator
  Teuchos::RCP<Structure> tmpstr = STR::TimIntCreate(*ioflags, *sdyn, *xparams, actdis, solver, contactsolver, output);
//  structure_ = STR::TimIntCreate(*ioflags, *sdyn, *xparams, actdis, solver, contactsolver, output);





  // original version
//  Teuchos::RCP<STR::TimInt> sti;
//  Teuchos::RCP<STR::TimIntImpl> stii;
//
//  sti = stii = STR::TimIntImplCreate(*ioflags, *sdyn, *xparams, actdis, solver, contactsolver, output);
//
//  Teuchos::RCP<STR::TimIntExpl> stie;
//  if (stii==Teuchos::null)
//  {
//    sti = stie = STR::TimIntExplCreate(*ioflags, *sdyn, *xparams, actdis, solver, contactsolver, output);
//  }






//  /* TODO: the following methods are still needed and will be reactivated as soon as possible!!! 13.03.2012 Georg Hammerl

  // still needed
//  // create auxiliar time integrator
//  Teuchos::RCP<STR::TimAda> sta = STR::TimAdaCreate(*ioflags, *sdyn, *xparams, *tap, sti);
//
//  if (sta!=Teuchos::null)
//  {
//    if (genprob.probtyp == prb_fsi or
//        genprob.probtyp == prb_fsi_lung or
//        genprob.probtyp == prb_gas_fsi or
//        genprob.probtyp == prb_biofilm_fsi or
//        genprob.probtyp == prb_thermo_fsi or
//        genprob.probtyp == prb_fluid_fluid_fsi)
//    {
//      dserror("no adaptive time integration with fsi");
//    }
//    structure_ = Teuchos::rcp(new StructureTimIntAda(sta, sti, ioflags, sdyn, xparams,
//                                                     actdis, solver, contactsolver, output));
//  }






  // no longer needed
//  else if (stie!=Teuchos::null)
//  {
//    if (genprob.probtyp == prb_fsi or
//        genprob.probtyp == prb_fsi_lung or
//        genprob.probtyp == prb_gas_fsi or
//        genprob.probtyp == prb_biofilm_fsi or
//        genprob.probtyp == prb_thermo_fsi or
//        genprob.probtyp == prb_fluid_fluid_fsi)
//    {
//      dserror("no explicit time integration with fsi");
//    }
//    structure_ = Teuchos::rcp(new StructureTimIntExpl(stie, ioflags, sdyn, xparams,
//                                                      actdis, solver, contactsolver, output));
//  }





//  else if (structure_ != Teuchos::null)
  if (tmpstr != Teuchos::null)
  {
    // no longer needed
//    Teuchos::RCP<Structure> tmpstr
//      = Teuchos::rcp(new StructureTimIntImpl(stii, ioflags, sdyn, xparams,
//                                             actdis, solver, contactsolver, output));

    if (genprob.probtyp == prb_fsi or
        genprob.probtyp == prb_fsi_lung or
        genprob.probtyp == prb_gas_fsi or
        genprob.probtyp == prb_biofilm_fsi or
        genprob.probtyp == prb_thermo_fsi or
        genprob.probtyp == prb_fsi_xfem or
        genprob.probtyp == prb_fluid_fluid_fsi)
    {
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

      Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);
      if ((actdis->Comm()).MyPID()==0) cout << "Using StructureNOXCorrectionWrapper()..." << endl;

      if (tmpstr->HaveConstraint())
      {
        if (coupling == fsi_iter_constr_monolithicstructuresplit or
            coupling == fsi_iter_constr_monolithicfluidsplit)
          structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
        else
          structure_ = rcp(new StructureNOXCorrectionWrapper(rcp(new StructureConstrMerged(tmpstr))));
      }
      else
      {
        if (coupling == fsi_iter_lung_monolithicstructuresplit or
            coupling == fsi_iter_lung_monolithicfluidsplit)
          structure_ = rcp(new StructureLung(rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        else
          structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
      }
    }
    else if (genprob.probtyp == prb_poroelast)
    {
      tmpstr = Teuchos::rcp(new StructureTimIntImplPoro(tmpstr));
      if (tmpstr->HaveConstraint())
      {
        structure_ = rcp(new StructureConstrMerged(tmpstr));
      }
      else
        structure_ = tmpstr;
    }
    else
    {
      structure_ = tmpstr;
    }
  }
  else
  {
    dserror("no proper time integration found");
  }





  // see you
  return;
}

Teuchos::RCP<LINALG::Solver> ADAPTER::StructureBaseAlgorithm::CreateLinearSolver(Teuchos::RCP<DRT::Discretization>& actdis)
{
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  solver = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                                    actdis->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  return solver;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> ADAPTER::StructureBaseAlgorithm::CreateContactMeshtyingSolver(Teuchos::RCP<DRT::Discretization>& actdis)
{
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  const Teuchos::ParameterList& mcparams     = DRT::Problem::Instance()->MeshtyingAndContactParams();
  INPAR::CONTACT::ApplicationType apptype    = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(mcparams,"APPLICATION");
  switch(DRT::INPUT::IntegralValue<int>(mcparams,"SYSTEM"))
  {
    case INPAR::CONTACT::system_spsimpler:
    {
      // meshtying/contact for structure
      if(apptype == INPAR::CONTACT::app_mortarmeshtying)
      {
        // plausibility check
        INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->MeshtyingSolverParams(),"AZPREC");
        if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
            prec != INPAR::SOLVER::azprec_TekoSIMPLE)
          dserror("Mortar/Contact with saddlepoint system only possible with SIMPLE preconditioner. Choose CheapSIMPLE or TekoSIMPLE in the MESHTYING SOLVER block in your dat file.");

        // build meshtying solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->MeshtyingSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      else if(apptype == INPAR::CONTACT::app_mortarcontact)
      {
        // plausibility check
        INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->ContactSolverParams(),"AZPREC");
        if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
            prec != INPAR::SOLVER::azprec_TekoSIMPLE)
          dserror("Mortar/Contact with saddlepoint system only possible with SIMPLE preconditioner. Choose CheapSIMPLE or TekoSIMPLE in the CONTACT SOLVER block in your dat file.");

        // build contact solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->ContactSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      else
      {
	// TODO: handle constraint solver case
        // plausibility check
        INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->ContactSolverParams(),"AZPREC");
        if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
            prec != INPAR::SOLVER::azprec_TekoSIMPLE)
          dserror("Mortar/Contact with saddlepoint system only possible with SIMPLE preconditioner. Choose CheapSIMPLE or TekoSIMPLE in the CONTACT SOLVER block in your dat file.");

        // build contact solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->ContactSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      actdis->ComputeNullSpaceIfNecessary(solver->Params());

      INPAR::CONTACT::ApplicationType apptype = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(mcparams,"APPLICATION");
      if     (apptype == INPAR::CONTACT::app_mortarcontact) solver->Params().set<bool>("CONTACT",true);
      else if(apptype==INPAR::CONTACT::app_mortarmeshtying) solver->Params().set<bool>("MESHTYING",true);
      else dserror("this cannot be: no saddlepoint problem for beamcontact or pure structure problem.");

      INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(mcparams,"STRATEGY");
      INPAR::CONTACT::SystemType      systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(mcparams,"SYSTEM");
      if (soltype==INPAR::CONTACT::solution_lagmult && systype!=INPAR::CONTACT::system_condensed)
      {
        solver->PutSolverParamsToSubParams("Inverse1", DRT::Problem::Instance()->StructSolverParams());
        solver->PutSolverParamsToSubParams("Inverse2", DRT::Problem::Instance()->ContactConstraintSolverParams());

        // note: the null space is definitely too long and wrong for the Lagrange multipliers
        // don't forget to call FixMLNullspace for "Inverse1"!
        actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1"));
      }
    }
    break;
    default:
    {
      if(apptype == INPAR::CONTACT::app_mortarmeshtying)
      {
        // build meshtying solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->MeshtyingSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      else if(apptype == INPAR::CONTACT::app_mortarcontact)
      {
        // build contact solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->ContactSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      else
      {
	// TODO: handle constraint solver case
        // build contact solver
        solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->ContactSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      }
      actdis->ComputeNullSpaceIfNecessary(solver->Params());
    }
  }

  return solver;
}

/*----------------------------------------------------------------------*/
#endif
