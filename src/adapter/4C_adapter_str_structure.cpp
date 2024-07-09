/*----------------------------------------------------------------------*/
/*! \file

\brief Deprecated structure field adapter
       (see ad_str_structure_new.H/.cpp for the new version)

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_str_structure.hpp"

#include "4C_adapter_str_constr_merged.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_adapter_str_fsi_timint_adaptive.hpp"
#include "4C_adapter_str_fsiwrapper_immersed.hpp"
#include "4C_adapter_str_lung.hpp"
#include "4C_adapter_str_redairway.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_timeloop.hpp"
#include "4C_adapter_str_timint_adaptive.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_poroelast.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_structure_timada_create.hpp"
#include "4C_structure_timint_create.hpp"
#include "4C_structure_timint_impl.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::StructureBaseAlgorithm::StructureBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& sdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  create_structure(prbdyn, sdyn, actdis);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithm::create_structure(const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& sdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  // major switch to different time integrators
  switch (Core::UTILS::IntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYP"))
  {
    case Inpar::Solid::dyna_statics:
    case Inpar::Solid::dyna_genalpha:
    case Inpar::Solid::dyna_onesteptheta:
    case Inpar::Solid::dyna_expleuler:
    case Inpar::Solid::dyna_centrdiff:
    case Inpar::Solid::dyna_ab2:
    case Inpar::Solid::dyna_euma:
    case Inpar::Solid::dyna_euimsto:
      create_tim_int(prbdyn, sdyn, actdis);  // <-- here is the show
      break;
    default:
      FOUR_C_THROW(
          "unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithm::create_tim_int(const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& sdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("Adapter::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // get the problem instance
  Global::Problem* problem = Global::Problem::instance();
  // what's the current problem type?
  Core::ProblemType probtype = problem->get_problem_type();

  // get mortar information
  std::vector<Core::Conditions::Condition*> mtcond(0);
  std::vector<Core::Conditions::Condition*> ccond(0);
  actdis->get_condition("Mortar", mtcond);
  actdis->get_condition("Contact", ccond);
  bool onlymeshtying = false;
  bool onlycontact = false;
  bool meshtyingandcontact = false;
  if (mtcond.size() != 0 and ccond.size() != 0) meshtyingandcontact = true;
  if (mtcond.size() != 0 and ccond.size() == 0) onlymeshtying = true;
  if (mtcond.size() == 0 and ccond.size() != 0) onlycontact = true;

  // Problem-types involving changing mesh or redistribution of mesh
  // for load balancing (like contact) during the simulation needs an additional step.
  // This is because the discretization read from the input file
  // do not match with the discr. at the current time step.
  // Here we read the discretization at the current time step from restart files
  if (not actdis->filled() || not actdis->have_dofs()) actdis->fill_complete();

  // get input parameter lists and copy them, because a few parameters are overwritten
  // const Teuchos::ParameterList& probtype
  //  = problem->ProblemTypeParams();
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->io_params()));
  Teuchos::RCP<Teuchos::ParameterList> tap =
      Teuchos::rcp(new Teuchos::ParameterList(sdyn.sublist("TIMEADAPTIVITY")));
  Teuchos::RCP<Teuchos::ParameterList> snox =
      Teuchos::rcp(new Teuchos::ParameterList(problem->structural_nox_params()));

  // show default parameters
  if ((actdis->get_comm()).MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::ParameterList& nox = xparams->sublist("NOX");
  nox = *snox;

  // Check if for chosen Rayleigh damping the regarding parameters are given explicitly in the .dat
  // file
  if (Core::UTILS::IntegralValue<Inpar::Solid::DampKind>(sdyn, "DAMPING") ==
      Inpar::Solid::damp_rayleigh)
  {
    if (sdyn.get<double>("K_DAMP") < 0.0)
    {
      FOUR_C_THROW("Rayleigh damping parameter K_DAMP not explicitly given.");
    }
    if (sdyn.get<double>("M_DAMP") < 0.0)
    {
      FOUR_C_THROW("Rayleigh damping parameter M_DAMP not explicitly given.");
    }
  }

  // create a solver
  Teuchos::RCP<Core::LinAlg::Solver> solver = create_linear_solver(actdis, sdyn);

  // create contact/meshtying solver only if contact/meshtying problem.
  Teuchos::RCP<Core::LinAlg::Solver> contactsolver = Teuchos::null;

  if (onlymeshtying or onlycontact or meshtyingandcontact)
    contactsolver = create_contact_meshtying_solver(actdis, sdyn);

  if (solver != Teuchos::null && (solver->params().isSublist("Belos Parameters")) &&
      solver->params().isSublist("ML Parameters")  // TODO what about MueLu?
      && Core::UTILS::IntegralValue<Inpar::Solid::StcScale>(sdyn, "STC_SCALING") !=
             Inpar::Solid::stc_none)
  {
    Teuchos::ParameterList& mllist = solver->params().sublist("ML Parameters");
    Teuchos::RCP<std::vector<double>> ns =
        mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace");

    const int size = actdis->dof_row_map()->NumMyElements();

    // extract the six nullspace vectors corresponding to the modes
    // trans x, trans y, trans z, rot x, rot y, rot z
    // Note: We assume 3d here!

    Teuchos::RCP<Epetra_Vector> nsv1 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), ns->data()));
    Teuchos::RCP<Epetra_Vector> nsv2 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), &ns->at(size)));
    Teuchos::RCP<Epetra_Vector> nsv3 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), &ns->at(2 * size)));
    Teuchos::RCP<Epetra_Vector> nsv4 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), &ns->at(3 * size)));
    Teuchos::RCP<Epetra_Vector> nsv5 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), &ns->at(4 * size)));
    Teuchos::RCP<Epetra_Vector> nsv6 =
        Teuchos::rcp(new Epetra_Vector(View, *(actdis->dof_row_map()), &ns->at(5 * size)));


    // prepare matrix for scaled thickness business of thin shell structures
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcinv =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*actdis->dof_row_map(), 81, true, true));

    stcinv->zero();
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // action for elements
    const std::string action = "calc_stc_matrix_inverse";
    p.set("action", action);
    p.set<int>(
        "stc_scaling", Core::UTILS::IntegralValue<Inpar::Solid::StcScale>(sdyn, "STC_SCALING"));
    p.set("stc_layer", 1);

    actdis->evaluate(p, stcinv, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

    stcinv->complete();

    for (int lay = 2; lay <= sdyn.get<int>("STC_LAYER"); ++lay)
    {
      Teuchos::ParameterList pe;

      p.set("stc_layer", lay);

      Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpstcmat =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*actdis->dof_row_map(), 81, true, true));
      tmpstcmat->zero();

      actdis->evaluate(p, tmpstcmat, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
      tmpstcmat->complete();

      stcinv = MLMultiply(*stcinv, *tmpstcmat, false, false, true);
    }

    Teuchos::RCP<Epetra_Vector> temp = Core::LinAlg::CreateVector(*(actdis->dof_row_map()), false);

    stcinv->multiply(false, *nsv1, *temp);
    nsv1->Update(1.0, *temp, 0.0);
    stcinv->multiply(false, *nsv2, *temp);
    nsv2->Update(1.0, *temp, 0.0);
    stcinv->multiply(false, *nsv3, *temp);
    nsv3->Update(1.0, *temp, 0.0);
    stcinv->multiply(false, *nsv4, *temp);
    nsv4->Update(1.0, *temp, 0.0);
    stcinv->multiply(false, *nsv5, *temp);
    nsv5->Update(1.0, *temp, 0.0);
    stcinv->multiply(false, *nsv6, *temp);
    nsv6->Update(1.0, *temp, 0.0);
  }

  // Checks in case of multi-scale simulations
  {
    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<Mat::PAR::Bundle> materials = problem->materials();
    for (const auto& [_, par] : materials->map())
    {
      if (par->type() == Core::Materials::m_struct_multiscale)
      {
        if (Core::UTILS::IntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYP") !=
            Inpar::Solid::dyna_genalpha)
          FOUR_C_THROW("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (Core::UTILS::IntegralValue<Inpar::Solid::MidAverageEnum>(
                     sdyn.sublist("GENALPHA"), "GENAVG") != Inpar::Solid::midavg_trlike)
          FOUR_C_THROW(
              "In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=TrLike");
        break;
      }
    }
  }

  // context for output and restart
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->writer();
  if (Core::UTILS::IntegralValue<int>(*ioflags, "OUTPUT_BIN"))
  {
    output->write_mesh(0, 0.0);
  }

  // create marching time integrator
  Teuchos::RCP<Solid::TimInt> tmpstr =
      Solid::TimIntCreate(prbdyn, *ioflags, sdyn, *xparams, actdis, solver, contactsolver, output);
  // initialize the time integrator
  tmpstr->init(prbdyn, sdyn, *xparams, actdis, solver);

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  /* Overwrite certain parameters in STRUCTURAL DYNAMIC/TIMADAPTIVITY by those from
   * FSI DYNAMIC/TIMEADAPTIVITY
   *
   * In case, that the structure field is part of an FSI simulation with time step
   * size adaptivity based on structure field error estimation, we have to provide
   * the following algorithmic control parameters:
   *
   * - ADAPTSTEPMAX
   * - STEPSIZEMAX
   * - STEPSIZEMIN
   * - SIZERATIOMAX
   * - SIZERATIOMIN
   * - SIZERATIOSCALE
   *
   * They are specified by the FSI algorithm since they have to be the same for
   * the structure and fluid field. Hence, we overwrite the corresponding
   * parameters in the structural parameter list in order to avoid redundant
   * parameter specification in the input file.
   *
   * Note: This is really ugly, but currently the only way to avoid that the user
   * has to specify these parameters twice in the input file.
   *
   * ToDO: Find something nicer here!
   *
   * \author mayr.mt \date 12/2013
   */
  // ---------------------------------------------------------------------------
  if (probtype == Core::ProblemType::fsi or probtype == Core::ProblemType::fsi_redmodels)
  {
    const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();
    const Teuchos::ParameterList& fsiada = fsidyn.sublist("TIMEADAPTIVITY");
    if (Core::UTILS::IntegralValue<bool>(fsiada, "TIMEADAPTON"))
    {
      // overrule time step size adaptivity control parameters
      if (tap->get<std::string>("KIND") != "NONE")
      {
        tap->set<int>("ADAPTSTEPMAX", fsiada.get<int>("ADAPTSTEPMAX"));
        tap->set<double>("STEPSIZEMAX", fsiada.get<double>("DTMAX"));
        tap->set<double>("STEPSIZEMIN", fsiada.get<double>("DTMIN"));
        tap->set<double>("SIZERATIOMAX", fsiada.get<double>("SIZERATIOMAX"));
        tap->set<double>("SIZERATIOMIN", fsiada.get<double>("SIZERATIOMIN"));
        tap->set<double>("SIZERATIOSCALE", fsiada.get<double>("SAFETYFACTOR"));

        if (actdis->get_comm().MyPID() == 0)
        {
          Core::IO::cout
              << "*** Due to FSI time step size adaptivity with structure based error estimation,\n"
                 "algorithmic control parameters in STRUCTURAL DYNAMIC/TIMEADAPTIVITY have been\n"
                 "overwritten by those from FSI DYNAMIC/TIMEADAPTIVITY."
              << Core::IO::endl
              << Core::IO::endl;
        }
      }
    }
  }
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  // create auxiliary time integrator, can be seen as a wrapper for tmpstr
  Teuchos::RCP<Solid::TimAda> sta =
      Solid::TimAdaCreate(*ioflags, prbdyn, sdyn, *xparams, *tap, tmpstr);

  if (sta != Teuchos::null and tmpstr != Teuchos::null)
  {
    switch (probtype)
    {
      case Core::ProblemType::structure:  // pure structural time adaptivity
      {
        structure_ = Teuchos::rcp(new StructureTimIntAda(sta, tmpstr));
        break;
      }
      case Core::ProblemType::fsi:  // structure based time adaptivity within an FSI simulation
      case Core::ProblemType::fsi_redmodels:
      {
        if ((actdis->get_comm()).MyPID() == 0)
          Core::IO::cout << "Using StructureNOXCorrectionWrapper()..." << Core::IO::endl;

        Teuchos::RCP<FSIStructureWrapper> fsiwrapperwithadaptivity =
            Teuchos::rcp(new StructureFSITimIntAda(
                sta, Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        // strTeuchos::rcp_dynamic_cast<StructureFSITimIntAda>(fsiwrapperwithadaptivity)->GetStrTimIntPtr();
        structure_ = fsiwrapperwithadaptivity;
        // structure_->GetStrTimIntPtr()-(prbdyn,sdyn,*xparams,actdis,solver);
        break;
      }
      default:
      {
        FOUR_C_THROW(
            "Adaptive time integration for the structure not implemented for desired problem "
            "type.");
        break;
      }
    }
  }
  else if (sta == Teuchos::null and tmpstr != Teuchos::null)
  {
    switch (probtype)
    {
      case Core::ProblemType::fsi:
      case Core::ProblemType::fsi_redmodels:
      case Core::ProblemType::fsi_lung:
      case Core::ProblemType::gas_fsi:
      case Core::ProblemType::ac_fsi:
      case Core::ProblemType::biofilm_fsi:
      case Core::ProblemType::thermo_fsi:
      {
        const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();
        const int coupling = Core::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");

        if ((actdis->get_comm()).MyPID() == 0)
          Core::IO::cout << "Using StructureNOXCorrectionWrapper()..." << Core::IO::endl;

        if (tmpstr->have_constraint())
        {
          if (coupling == fsi_iter_constr_monolithicstructuresplit or
              coupling == fsi_iter_constr_monolithicfluidsplit)
            structure_ = Teuchos::rcp(
                new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
          else
            structure_ = Teuchos::rcp(
                new StructureConstrMerged(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        }
        else
        {
          if (coupling == fsi_iter_lung_monolithicstructuresplit or
              coupling == fsi_iter_lung_monolithicfluidsplit)
            structure_ = Teuchos::rcp(
                new StructureLung(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
          else
            structure_ = Teuchos::rcp(
                new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(tmpstr))));
        }
      }
      break;
      case Core::ProblemType::immersed_fsi:
      {
        structure_ = Teuchos::rcp(new FSIStructureWrapperImmersed(tmpstr));
      }
      break;
      case Core::ProblemType::ssi:
      case Core::ProblemType::ssti:
      {
        structure_ = Teuchos::rcp(new SSIStructureWrapper(tmpstr));
      }
      break;
      case Core::ProblemType::redairways_tissue:
      {
        structure_ = Teuchos::rcp(new StructureRedAirway(tmpstr));
      }
      break;
      case Core::ProblemType::poroelast:
      case Core::ProblemType::poroscatra:
      case Core::ProblemType::fpsi:
      case Core::ProblemType::fps3i:
      case Core::ProblemType::fpsi_xfem:
      case Core::ProblemType::fsi_xfem:
      {
        const Teuchos::ParameterList& porodyn = problem->poroelast_dynamic_params();
        const Inpar::PoroElast::SolutionSchemeOverFields coupling =
            Core::UTILS::IntegralValue<Inpar::PoroElast::SolutionSchemeOverFields>(
                porodyn, "COUPALGO");
        if (tmpstr->have_constraint())
        {
          if (coupling == Inpar::PoroElast::Monolithic_structuresplit or
              coupling == Inpar::PoroElast::Monolithic_fluidsplit or
              coupling == Inpar::PoroElast::Monolithic_nopenetrationsplit)
            structure_ = Teuchos::rcp(new FPSIStructureWrapper(tmpstr));
          else
            structure_ = Teuchos::rcp(new StructureConstrMerged(tmpstr));
        }
        else
        {
          structure_ = Teuchos::rcp(new FPSIStructureWrapper(tmpstr));
        }
      }
      break;
      case Core::ProblemType::struct_ale:
      {
        structure_ = Teuchos::rcp(new FSIStructureWrapper(tmpstr));
      }
      break;
      default:
      {
        /// wrap time loop for pure structure problems
        structure_ = (Teuchos::rcp(new StructureTimeLoop(tmpstr)));
      }
      break;
    }
  }
  else
  {
    FOUR_C_THROW("no proper time integration found");
  }
  // see you
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> Adapter::StructureBaseAlgorithm::create_linear_solver(
    Teuchos::RCP<Core::FE::Discretization>& actdis, const Teuchos::ParameterList& sdyn)
{
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::null;

  // get the solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  solver = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::instance()->solver_params(linsolvernumber),
          actdis->get_comm(), Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY")));

  actdis->compute_null_space_if_necessary(solver->params());

  return solver;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> Adapter::StructureBaseAlgorithm::create_contact_meshtying_solver(
    Teuchos::RCP<Core::FE::Discretization>& actdis, const Teuchos::ParameterList& sdyn)
{
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::null;

  // Get mortar information: contact or meshtying or both?
  bool onlymeshtying = false;
  bool onlycontact = false;
  bool meshtyingandcontact = false;
  {
    std::vector<Core::Conditions::Condition*> mtcond(0);
    std::vector<Core::Conditions::Condition*> ccond(0);
    actdis->get_condition("Mortar", mtcond);
    actdis->get_condition("Contact", ccond);
    if (mtcond.size() != 0 and ccond.size() != 0) meshtyingandcontact = true;
    if (mtcond.size() != 0 and ccond.size() == 0) onlymeshtying = true;
    if (mtcond.size() == 0 and ccond.size() != 0) onlycontact = true;
  }
  const Teuchos::ParameterList& mcparams = Global::Problem::instance()->contact_dynamic_params();

  // Get the solver number used for meshtying/contact problems
  const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
  // check if the meshtying/contact solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "No linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in "
        "CONTACT DYNAMIC to a valid number!");

  // Distinguish the system type, i.e. condensed vs. saddle-point
  switch (Core::UTILS::IntegralValue<int>(mcparams, "SYSTEM"))
  {
    case Inpar::CONTACT::system_saddlepoint:
    {
      /* Plausibility check
       *
       * Solver can be either a direct solver (UMFPACK, Superlu) or an iterative solver (Belos).
       */
      const auto sol = Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(
          Global::Problem::instance()->solver_params(linsolvernumber), "SOLVER");
      const auto prec = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
          Global::Problem::instance()->solver_params(linsolvernumber), "AZPREC");
      if (sol != Core::LinearSolver::SolverType::umfpack &&
          sol != Core::LinearSolver::SolverType::superlu)
      {
        // if an iterative solver is chosen we need a block preconditioner
        if (prec != Core::LinearSolver::PreconditionerType::cheap_simple &&
            prec != Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp)
          FOUR_C_THROW(
              "You have chosen an iterative linear solver. For mortar meshtying/contact problems "
              "in saddle-point formulation, a block preconditioner is required. Choose an "
              "appropriate block preconditioner such as CheapSIMPLE or MueLu_contactSP "
              "(if MueLu is available) in the SOLVER %i block in your input file.",
              linsolvernumber);
      }

      // build meshtying/contact solver
      solver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::instance()->solver_params(linsolvernumber),
              actdis->get_comm(), Global::Problem::instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::instance()->io_params(), "VERBOSITY")));

      actdis->compute_null_space_if_necessary(solver->params());

      // feed the solver object with additional information
      if (onlycontact or meshtyingandcontact)
        solver->params().set<bool>("CONTACT", true);
      else if (onlymeshtying)
        solver->params().set<bool>("MESHTYING", true);
      else
        FOUR_C_THROW(
            "Saddle-point formulations are only supported for solid CONTACT or MESHTYING problems. "
            "Problems like beamcontact or pure structure problem w/o contact do not support a "
            "saddle-point formulation.");

      Inpar::CONTACT::SolvingStrategy soltype =
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(mcparams, "STRATEGY");
      if (soltype == Inpar::CONTACT::solution_lagmult)
      {
        // get the solver number used for structural problems
        const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
        // check if the structural solver has a valid solver number
        if (linsolvernumber == (-1))
          FOUR_C_THROW(
              "No linear solver defined for structural field. Please set LINEAR_SOLVER in "
              "STRUCTURAL DYNAMIC to a valid number!");

        // provide null space information
        if (prec == Core::LinearSolver::PreconditionerType::cheap_simple)
        {
          actdis->compute_null_space_if_necessary(
              solver->params()
                  .sublist("CheapSIMPLE Parameters")
                  .sublist("Inverse1"));  // Inverse2 is created within blockpreconditioners.cpp
        }
        else if (prec == Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp)
        { /* do nothing here */
        }
      }
    }
    break;
    default:
    {
      // build meshtying solver
      solver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::instance()->solver_params(linsolvernumber),
              actdis->get_comm(), Global::Problem::instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::instance()->io_params(), "VERBOSITY")));
      actdis->compute_null_space_if_necessary(solver->params());
    }
    break;
  }

  return solver;
}

FOUR_C_NAMESPACE_CLOSE
