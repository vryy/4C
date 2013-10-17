/*----------------------------------------------------------------------*/
/*!
\file strtimint.cpp
\brief Time integration for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "../drt_io/io_pstream.H"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_TimeMonitor.hpp"

#include "strtimint_mstep.H"
#include "strtimint.H"
#include "stru_resulttest.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/micromaterial.H"

#include "../drt_lib/drt_locsys.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_potential/drt_potential_manager.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_contact/contact_analytical.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_abstract_strategy.H" // for feeding contactsolver with maps
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_patspec/patspec.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_plastic_ssn/plastic_ssn_manager.H"
#include "../drt_stru_multi/microstatic.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_solver.H"

#include "../drt_so3/so_sh8p8.H"
#include "../drt_so3/so3_ssn_plast_eletypes.H"
#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*/
/* print tea time logo */
void STR::TimInt::Logo()
{
 IO::cout << "Welcome to Structural Time Integration " <<IO::endl;
 IO::cout << "     __o__                          __o__" <<IO::endl;
 IO::cout << "__  /-----\\__                  __  /-----\\__" <<IO::endl;
 IO::cout << "\\ \\/       \\ \\    |       \\    \\ \\/       \\ \\" <<IO::endl;
 IO::cout << " \\ |  tea  | |    |-------->    \\ |  tea  | |" <<IO::endl;
 IO::cout << "  \\|       |_/    |       /      \\|       |_/" <<IO::endl;
 IO::cout << "    \\_____/   ._                   \\_____/   ._ _|_ /|" <<IO::endl;
 IO::cout << "              | |                            | | |   |" <<IO::endl;
 IO::cout <<IO::endl;
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimInt::TimInt
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: discret_(actdis),
  myrank_(actdis->Comm().MyPID()),
  dofrowmap_(actdis->Filled() ? actdis->DofRowMap() : NULL),
  solver_(solver),
  contactsolver_(contactsolver),
  solveradapttol_(DRT::INPUT::IntegralValue<int>(sdynparams,"ADAPTCONV")==1),
  solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
  dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor())),
  divcontype_(DRT::INPUT::IntegralValue<INPAR::STR::DivContAct>(sdynparams,"DIVERCONT")),
  output_(output),
  printlogo_(true),  // DON'T EVEN DARE TO SET THIS TO FALSE
  printscreen_(ioparams.get<int>("STDOUTEVRY")),
  errfile_(xparams.get<FILE*>("err file")),
  printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
  printiter_(true),  // ADD INPUT PARAMETER
  writerestartevery_(sdynparams.get<int>("RESTARTEVRY")),
  writereducedrestart_(xparams.get<int>("REDUCED_OUTPUT")),
  writestate_((bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_DISP")),
  writevelacc_((bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_VEL_ACC")),
  writeresultsevery_(sdynparams.get<int>("RESULTSEVRY")),
  writestress_(DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_STRESS")),
  writecouplstress_(DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_COUPLING_STRESS")),
  writestrain_(DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_STRAIN")),
  writeplstrain_(DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_PLASTIC_STRAIN")),
  writeenergyevery_(sdynparams.get<int>("RESEVRYERGY")),
  writesurfactant_((bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_SURFACTANT")),
  energyfile_(NULL),
  damping_(DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(sdynparams,"DAMPING")),
  dampk_(sdynparams.get<double>("K_DAMP")),
  dampm_(sdynparams.get<double>("M_DAMP")),
  conman_(Teuchos::null),
  consolv_(Teuchos::null),
  surfstressman_(Teuchos::null),
  potman_(Teuchos::null),
  cmtman_(Teuchos::null),
  beamcman_(Teuchos::null),
  statmechman_(Teuchos::null),
  locsysman_(Teuchos::null),
  pressure_(Teuchos::null),
  time_(Teuchos::null),
  timen_(0.0),
  dt_(Teuchos::null),
  timemax_(sdynparams.get<double>("MAXTIME")),
  stepmax_(sdynparams.get<int>("NUMSTEP")),
  step_(0),
  stepn_(0),
  firstoutputofrun_(true),
  lumpmass_(DRT::INPUT::IntegralValue<int>(sdynparams,"LUMPMASS")==1),
  young_temp_(DRT::INPUT::IntegralValue<int>(sdynparams,"YOUNG_IS_TEMP_DEPENDENT")==1),
  zeros_(Teuchos::null),
  dis_(Teuchos::null),
  dism_(Teuchos::null),
  vel_(Teuchos::null),
  acc_(Teuchos::null),
  disn_(Teuchos::null),
  dismatn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  fifc_(Teuchos::null),
  stiff_(Teuchos::null),
  mass_(Teuchos::null),
  damp_(Teuchos::null),
  timer_(Teuchos::rcp(new Epetra_Time(actdis->Comm()))),
  dtsolve_(0.0),
  dtele_(0.0),
  dtcmt_(0.0),
  couplstate_(0),
  pslist_(Teuchos::null),
  strgrdisp_(Teuchos::null)
{
  // welcome user
  if ( (printlogo_) and (myrank_ == 0) )
  {
    Logo();
  }

  // check whether discretisation has been completed
  if (not discret_->Filled() || not actdis->HaveDofs())
  {
    dserror("Discretisation is not complete or has no dofs!");
  }

  // connect degrees of freedom for periodic boundary conditions
  // (i.e. for multi-patch nurbs computations)
  {
    PeriodicBoundaryConditions pbc(discret_);

    if (pbc.HasPBC())
    {
      pbc.UpdateDofsForPeriodicBoundaryConditions();

      discret_->ComputeNullSpaceIfNecessary(solver->Params(),true);

      dofrowmap_ = discret_->DofRowMap(0);
    }
  }

  // time state
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, 0.0));  // HERE SHOULD BE SOMETHING LIKE (sdynparams.get<double>("TIMEINIT"))
  dt_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, sdynparams.get<double>("TIMESTEP")));
  step_ = 0;
  timen_ = (*time_)[0] + (*dt_)[0];  // set target time to initial time plus step size
  stepn_ = step_ + 1;

  // output file for energy
  if ( (writeenergyevery_ != 0) and (myrank_ == 0) )
    AttachEnergyFile();

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap_, true);

  // Map containing Dirichlet DOFs
  {
    Teuchos::ParameterList p;
    p.set("total time", timen_);
    discret_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  // displacements D_{n}
  dis_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));

  // displacements D_{n+1} at t_{n+1}
  disn_ = LINALG::CreateVector(*dofrowmap_, true);

  if (DRT::Problem::Instance()->ProblemType() == prb_struct_ale and
      (DRT::Problem::Instance()->ContactDynamicParams()).get<double>("WEARCOEFF")>0.0)
  {
    // material displacements Dm_{n+1} at t_{n+1}
    dismatn_ = LINALG::CreateVector(*dofrowmap_,true);

    // material_displacements D_{n}
    dism_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  }

  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*dofrowmap_, true);
  // accelerations A_{n+1} at t_{n+1}
  accn_ = LINALG::CreateVector(*dofrowmap_, true);
  // create empty interface force vector
  fifc_ = LINALG::CreateVector(*dofrowmap_, true);

  // create empty matrices
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));
  mass_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));
  if (damping_ != INPAR::STR::damp_none)
  {
    if (!HaveNonlinearMass())
    {
      damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));
    }
    else
    {
      //Since our element evaluate routine is only designed for two input matrices
      //(stiffness and damping or stiffness and mass) its not possible, to have nonlinear
      //inertia forces AND material damping.
      dserror("So far its not possible to model nonlinear inertia forces and damping!");
    }
  }


  // set initial fields
  SetInitialFields();

  // initialize constraint manager
  conman_ = Teuchos::rcp(new UTILS::ConstrManager(discret_,
                                                  (*dis_)(0),
                                                  sdynparams));
  // initialize constraint solver iff constraints are defined
  if (conman_->HaveConstraint())
  {
    consolv_ = Teuchos::rcp(new UTILS::ConstraintSolver(discret_,
                                                        *solver_,
                                                        dbcmaps_,
                                                        sdynparams));
  }

  // check for beam contact
  {
    // If beam contact (no statistical mechanics) is chosen in the input file, then a
    // corresponding manager object stored via #beamcman_ is created and all relevant
    // stuff is initialized. Else, #beamcman_ remains a Teuchos::null pointer.
    PrepareBeamContact(sdynparams);
  }

  // check for mortar contact or meshtying
  {
    // If mortar contact or meshtying is chosen in the input file, then a
    // corresponding manager object stored via #cmtman_ is created and all relevant
    // stuff is initialized. Else, #cmtman_ remains a Teuchos::null pointer.
    PrepareContactMeshtying(sdynparams);
  }

  // check for elements using a semi-smooth Newton method for plasticity
  {
    // If at least one such element exists, then a
    // corresponding manager object stored via #ssnplastman_ is created and all relevant
    // stuff is initialized. Else, #ssnplastman_ remains a Teuchos::null pointer.
    PrepareSemiSmoothPlasticity();
  }

  // check for statistical mechanics
  {
    PrepareStatMech();
  }
  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  // (this is necessary to BOTH constraints and contact / meshtying)
  dofrowmap_ = discret_->DofRowMap();

  // Initialize SurfStressManager for handling surface stress conditions due to interfacial phenomena
  surfstressman_ = Teuchos::rcp(new UTILS::SurfStressManager(discret_,
                                                             sdynparams,
                                                             DRT::Problem::Instance()->OutputControlFile()->FileName()));

  // Check for potential conditions
  {
    std::vector<DRT::Condition*> potentialcond(0);
    discret_->GetCondition("Potential", potentialcond);
    if (potentialcond.size())
    {
      potman_ = Teuchos::rcp(new POTENTIAL::PotentialManager(Discretization(), *discret_));
      stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false, LINALG::SparseMatrix::FE_MATRIX));
    }
  }

  // check whether we have locsys BCs and create LocSysManager if so
  // after checking
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      locsysman_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_));
    }
  }

  // check if we have elements which use a continuous displacement and pressure
  // field
  {
    int locnumsosh8p8 = 0;
    // Loop through all elements on processor
    for (int i=0; i<discret_->NumMyColElements(); ++i)
    {
      // get the actual element

      if (discret_->lColElement(i)->ElementType() == DRT::ELEMENTS::So_sh8p8Type::Instance())
        locnumsosh8p8 += 1;

    }
    // Was at least one SoSh8P8 found on one processor?
    int glonumsosh8p8 = 0;
    discret_->Comm().MaxAll(&locnumsosh8p8, &glonumsosh8p8, 1);
    // Yes, it was. Go ahead for all processors (even if they do not carry any SoSh8P8 elements)
    if (glonumsosh8p8 > 0)
    {
      pressure_ = Teuchos::rcp(new LINALG::MapExtractor());
      const int ndim = 3;
      FLD::UTILS::SetupFluidSplit(*discret_, ndim, *pressure_);
    }
  }

  // check if we have elements with a micro-material
  havemicromat_ = false;
  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele = discret_->lColElement(i);
    RCP<MAT::Material> mat = actele->Material();
    if (mat != Teuchos::null && mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      havemicromat_ = true;
      break;
    }
  }

  // check for patient specific needs
  const Teuchos::ParameterList& patspec  = DRT::Problem::Instance()->PatSpecParams();
  if (DRT::INPUT::IntegralValue<int>(patspec,"PATSPEC"))
  {
	pslist_ = Teuchos::rcp(new Teuchos::ParameterList());
	pslist_->set("haveembedtissue", false);
    // check if patspeccond are already initialized
    // this is of relevance for Montecarlo Simulation
    std::vector<DRT::Condition*> pscond;
    discret_->GetCondition("PatientSpecificData", pscond);
    if (!pscond.size())
    {
      //initialize patient specific parameters and conditions
      PATSPEC::PatientSpecificGeometry(discret_, pslist_);

      // fix pointer to dofrowmap_
      dofrowmap_ = discret_->DofRowMap();
    }
  }

  porositysplitter_ = POROELAST::UTILS::BuildPoroSplitter(discret_);

  // stay with us
  return;
}

/*----------------------------------------------------------------------*/
/* Set intitial fields in structure (e.g. initial velocities */
void STR::TimInt::SetInitialFields()
{
  //***************************************************
  // Data that needs to be handed into discretization:
  // - std::string field: name of initial field to be set
  // - std::vector<int> localdofs: local dof ids affected
  //***************************************************

  // set initial velocity field if existing
  const std::string field = "Velocity";
  std::vector<int> localdofs;
  localdofs.push_back(0);
  localdofs.push_back(1);
  localdofs.push_back(2);
  discret_->EvaluateInitialField(field,(*vel_)(0),localdofs);

  // set initial porosity field if existing
  const std::string porosityfield = "Porosity";
  std::vector<int> porositylocaldofs;
  porositylocaldofs.push_back(DRT::Problem::Instance()->NDim());
  discret_->EvaluateInitialField(porosityfield,(*dis_)(0),porositylocaldofs);

  return;
}

/*----------------------------------------------------------------------*/
/* Check for beam contact and do preparations */
void STR::TimInt::PrepareBeamContact(const Teuchos::ParameterList& sdynparams)
{
  // some parameters
  const Teuchos::ParameterList&   scontact = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::CONTACT::ApplicationType apptype  = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION");

  // only continue if beam contact unmistakably chosen in input file
  if (apptype == INPAR::CONTACT::app_beamcontact)
  {
    // store integration parameter alphaf into beamcman_ as well
    // (for all cases except OST, GenAlpha and GEMM this is zero)
    // (note that we want to hand in theta in the OST case, which
    // is defined just the other way round as alphaf in GenAlpha schemes.
    // Thus, we have to hand in 1.0-theta for OST!!!)
    double alphaf = 0.0;
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_genalpha)
      alphaf = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_gemm)
      alphaf = sdynparams.sublist("GEMM").get<double>("ALPHA_F");
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_onesteptheta)
      alphaf = 1.0 - sdynparams.sublist("ONESTEPTHETA").get<double>("THETA");

    // create beam contact manager
    beamcman_ = Teuchos::rcp(new CONTACT::Beam3cmanager(*discret_,alphaf));
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Check for contact or meshtying and do preparations */
void STR::TimInt::PrepareContactMeshtying(const Teuchos::ParameterList& sdynparams)
{
  // some parameters
  const Teuchos::ParameterList&   smortar  = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList&   scontact = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::MORTAR::ShapeFcn         shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar,"SHAPEFCN");
  INPAR::CONTACT::ApplicationType apptype  = DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION");
  INPAR::CONTACT::SolvingStrategy soltype  = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(scontact,"STRATEGY");
  bool semismooth = DRT::INPUT::IntegralValue<int>(scontact,"SEMI_SMOOTH_NEWTON");

  // check mortar contact or meshtying conditions
  std::vector<DRT::Condition*> mortarconditions(0);
  discret_->GetCondition("Mortar",mortarconditions);

  // only continue if mortar contact / meshtying unmistakably chosen in input file
  if (apptype == INPAR::CONTACT::app_mortarcontact || apptype == INPAR::CONTACT::app_mortarmeshtying)
  {
    // double-check for contact conditions
    if ((int)mortarconditions.size()==0)
      dserror("ERROR: No contact conditions provided, check your input file!");

    // store integration parameter alphaf into cmtman_ as well
    // (for all cases except OST, GenAlpha and GEMM this is zero)
    // (note that we want to hand in theta in the OST case, which
    // is defined just the other way round as alphaf in GenAlpha schemes.
    // Thus, we have to hand in 1-theta for OST!!!)
    double alphaf = 0.0;
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_genalpha)
      alphaf = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_gemm)
      alphaf = sdynparams.sublist("GEMM").get<double>("ALPHA_F");
    if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") == INPAR::STR::dyna_onesteptheta)
      alphaf = 1.0 - sdynparams.sublist("ONESTEPTHETA").get<double>("THETA");

    // decide whether this is meshtying or contact and create manager
    if (apptype == INPAR::CONTACT::app_mortarmeshtying)
      cmtman_ = Teuchos::rcp(new CONTACT::MtManager(*discret_,alphaf));
    else if (apptype == INPAR::CONTACT::app_mortarcontact)
      cmtman_ = Teuchos::rcp(new CONTACT::CoManager(*discret_,alphaf));

    // store DBC status in contact nodes
    cmtman_->GetStrategy().StoreDirichletStatus(dbcmaps_);

    // create old style dirichtoggle vector (supposed to go away)
    dirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->FullMap())));
    RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->CondMap())));
    temp->PutScalar(1.0);
    LINALG::Export(*temp,*dirichtoggle_);

    // contact and constraints together not yet implemented
    if (conman_->HaveConstraint())
      dserror("ERROR: Constraints and contact cannot be treated at the same time yet");

    // print messages for multifield problems (e.g FSI)
    const PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    const std::string probname = DRT::Problem::Instance()->ProblemName();
    if (probtype != prb_structure && !myrank_)
    {
      // warnings
#ifdef CONTACTPSEUDO2D
      std::cout << RED << "WARNING: The flag CONTACTPSEUDO2D is switched on. If this "
           << "is a real 3D problem, switch it off!" << END_COLOR << std::endl;
#else
      std::cout << RED << "WARNING: The flag CONTACTPSEUDO2D is switched off. If this "
           << "is a 2D problem modeled pseudo-3D, switch it on!" << END_COLOR << std::endl;
#endif // #ifdef CONTACTPSEUDO2D

//      if (probtype!=prb_tsi)
//        std::cout << RED << "WARNING: Contact and Meshtying are still experimental "
//             << "for the chosen problem type \"" << probname << "\"!\n" << END_COLOR << std::endl;

      // errors
      if (probtype!=prb_tsi and probtype!=prb_struct_ale)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && (!semismooth || shapefcn != INPAR::MORTAR::shape_dual))
          dserror("ERROR: Multifield problems with LM strategy for meshtying/contact only for dual+semismooth case!");
        if (soltype == INPAR::CONTACT::solution_auglag)
          dserror("ERROR: Multifield problems with AL strategy for meshtying/contact not yet implemented");
      }
    }

    // set zero displacment state
    cmtman_->GetStrategy().SetState("displacement",zeros_);

    // visualization of initial configuration
#ifdef MORTARGMSH3
    cmtman_->GetStrategy().VisualizeGmsh(0,0);
#endif // #ifdef MORTARGMSH3

    // initialization of contact or meshtying
    {
      // FOR MESHTYING (ONLY ONCE), NO FUNCTIONALITY FOR CONTACT CASES
      // (1) Do mortar coupling in reference configuration and
      // perform mesh initialization for rotational invariance
      cmtman_->GetStrategy().MortarCoupling(zeros_);
      cmtman_->GetStrategy().MeshInitialization();

      // FOR PENALTY CONTACT (ONLY ONCE), NO FUNCTIONALITY FOR OTHER CASES
      // (1) Explicitly store gap-scaling factor kappa
      cmtman_->GetStrategy().SaveReferenceState(zeros_);
    }

    //**********************************************************************
    // prepare solvers for contact/meshtying problem
    //**********************************************************************
    {
      // only plausability check, that a contact solver is available
      if (contactsolver_ == Teuchos::null)
        dserror("no contact solver in STR::TimInt::PrepareContactMeshtying? cannot be");
    }

    //**********************************************************************
    // feed solver/preconditioner with additional information about the contact/meshtying problem
    //**********************************************************************
    {
#if 0 // do we need this? feed solvers with latest information. why not using Aztec parameters?
      if (contactsolver_->Params().isSublist("MueLu (Contact) Parameters"))
      {
        Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("MueLu (Contact) Parameters");
        RCP<Epetra_Map> masterDofMap;
        RCP<Epetra_Map> slaveDofMap;
        RCP<Epetra_Map> innerDofMap;
        RCP<Epetra_Map> activeDofMap;
        // transform cmtman_ to CoAbstractStrategy object, since this code is only meant to work with contact/meshtying)
        Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtman_->GetStrategy());
        //Teuchos::RCP<CONTACT::CoAbstractStrategy> cstrat = Teuchos::rcp_dynamic_cast<CONTACT::CoAbstractStrategy>(strat);
        //if(cstrat != Teuchos::null) { // dserror("STR::TimInt::PrepareContactMeshtying: dynamic cast to CONTACT::CoAbstractStrategy failed. Are you running a contact/meshtying problem?");
        strat->CollectMapsForPreconditioner(masterDofMap, slaveDofMap, innerDofMap, activeDofMap);
        mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap",masterDofMap);
        mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap",slaveDofMap);
        mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::InnerDofMap",innerDofMap);
        mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap",activeDofMap);
        //} else std::cout << "STR::TimInt::PrepareContactMeshtying: dynamic cast to CONTACT::CoAbstractStrategy failed. Are you running a contact/meshtying problem? strtimint.cpp line 550" << std::endl;

        //std::cout << contactsolver_->Params() << std::endl;
      }

      // TODO fix me
      if (contactsolver_->Params().isSublist("MueLu (Contact2) Parameters"))
      {
        Teuchos::ParameterList& mueluParams = contactsolver_->Params().sublist("MueLu (Contact2) Parameters");
        RCP<Epetra_Map> masterDofMap;
        RCP<Epetra_Map> slaveDofMap;
        RCP<Epetra_Map> innerDofMap;
        RCP<Epetra_Map> activeDofMap;
        // transform cmtman_ to CoAbstractStrategy object, since this code is only meant to work with contact/meshtying)
        Teuchos::RCP<MORTAR::StrategyBase> strat = Teuchos::rcpFromRef(cmtman_->GetStrategy());
        //Teuchos::RCP<CONTACT::CoAbstractStrategy> cstrat = Teuchos::rcp_dynamic_cast<CONTACT::CoAbstractStrategy>(strat);
        //if(cstrat != Teuchos::null) { //dserror("STR::TimInt::PrepareContactMeshtying: dynamic cast to CONTACT::CoAbstractStrategy failed. Are you running a contact/meshtying problem?");
          strat->CollectMapsForPreconditioner(masterDofMap, slaveDofMap, innerDofMap, activeDofMap);
          mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap",masterDofMap);
          mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap",slaveDofMap);
          mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::InnerDofMap",innerDofMap);
          mueluParams.set<RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap",activeDofMap);
        //} else std::cout << "STR::TimInt::PrepareContactMeshtying: dynamic cast to CONTACT::CoAbstractStrategy failed. Are you running a contact/meshtying problem? strtimint.cpp line 550" << std::endl;
        //std::cout << contactsolver_->Params() << std::endl;
      }
#endif

    }

    // output of strategy type to screen
    {
      // output
      if (!myrank_)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_standard)
          std::cout << "===== Standard Lagrange multiplier strategy ====================\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
          std::cout << "===== Dual Lagrange multiplier strategy ========================\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
          std::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_standard)
          std::cout << "===== Standard Penalty strategy ================================\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_dual)
          std::cout << "===== Dual Penalty strategy ====================================\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_auglag && shapefcn == INPAR::MORTAR::shape_standard)
          std::cout << "===== Standard Augmented Lagrange strategy =====================\n" << std::endl;
        else if (soltype == INPAR::CONTACT::solution_auglag && shapefcn == INPAR::MORTAR::shape_dual)
          std::cout << "===== Dual Augmented Lagrange strategy =========================\n" << std::endl;
      }
    }
  }

  return;
}
/*----------------------------------------------------------------------*/
/* calculate stresses and strains on micro-scale */
void STR::TimInt::PrepareSemiSmoothPlasticity()
{
  int HavePlasticity_local=0;
  int HavePlasticity_global=0;
  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele = discret_->lColElement(i);
    if (   actele->ElementType() == DRT::ELEMENTS::So_hex8PlastType::Instance()
        || actele->ElementType() == DRT::ELEMENTS::So_tet4PlastType::Instance()
        || actele->ElementType() == DRT::ELEMENTS::So_hex8fbarPlastType::Instance()
       )
    {
      HavePlasticity_local=1;
      break;
    }
  }
  discret_->Comm().MaxAll(&HavePlasticity_local,&HavePlasticity_global,1);
  if (HavePlasticity_global)
    plastman_=Teuchos::rcp(new UTILS::PlastSsnManager(discret_));
  return;
}

/*----------------------------------------------------------------------*/
/* Prepare contact for new time step */
void STR::TimInt::PrepareStepContact()
{
  // just do something here if contact is present
  if (HaveContactMeshtying())
  {
    // set inttime_ to zero
    cmtman_->GetStrategy().Inttime_init();

    // dynamic parallel redistribution of interfaces
    cmtman_->GetStrategy().RedistributeContact((*dis_)(0));

    // evaluation of reference state for friction (only at t=0)
    cmtman_->GetStrategy().EvaluateReferenceState(step_,disn_);
  }

  // bye bye
  return;
}

/*----------------------------------------------------------------------*/
/* Check for contact or meshtying and do preparations */
void STR::TimInt::PrepareStatMech()
{
  // some parameters
  const Teuchos::ParameterList&   statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  INPAR::STATMECH::ThermalBathType tbtype  = DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams,"THERMALBATH");

  if(tbtype != INPAR::STATMECH::thermalbath_none)
  {
    statmechman_ = Teuchos::rcp(new STATMECH::StatMechManager(discret_));

    dirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
    RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->CondMap())));
    temp->PutScalar(1.0);
    LINALG::Export(*temp,*dirichtoggle_);

    // output
    if (!discret_->Comm().MyPID())
    {
      switch(tbtype)
      {
        case INPAR::STATMECH::thermalbath_uniform:
          std::cout << "========= Statistical Mechanics: uniform thermal bath ==========\n" << std::endl;
          break;
        case INPAR::STATMECH::thermalbath_shearflow:
          std::cout << "======== Statistical Mechanics: thermal bath, shearflow ========\n" << std::endl;
          break;
        default: dserror("Undefined thermalbath type!");
        break;
      }
      Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
      if(!ioparams.get<int>("STDOUTEVRY",0))
        std::cout<<"STDOUT SUPPRESSED!"<<std::endl;
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/* EvaluateReferenceState for frictional contact */
void STR::TimInt::EvaluateReferenceState()
{
  // only do something if contact is present
  if (HaveContactMeshtying())
    cmtman_->GetStrategy().EvaluateReferenceState(step_,disn_);
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state
 * and identify consistent accelerations */
void STR::TimInt::DetermineMassDampConsistAccel()
{
  // temporary force vectors in this routine
  Teuchos::RCP<Epetra_Vector> fext
    = LINALG::CreateVector(*dofrowmap_, true); // external force
  Teuchos::RCP<Epetra_Vector> fint
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  // initialise matrices
  stiff_->Zero();
  mass_->Zero();

  // overwrite initial state vectors with DirichletBCs
  ApplyDirichletBC((*time_)[0], (*dis_)(0), (*vel_)(0), (*acc_)(0), false);

  // get external force
  ApplyForceExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext, stiff_);

  // get initial internal force and stiffness and mass
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    if(lumpmass_ == false)
      p.set("action", "calc_struct_nlnstiffmass");
    // lumping the mass matrix
    else
      p.set("action", "calc_struct_nlnstifflmass");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    p.set<int>("young_temp", young_temp_);
    if (pressure_ != Teuchos::null) p.set("volume", 0.0);
    // set vector values needed by elements
    discret_->ClearState();
    // extended SetState(0,...) in case of multiple dofsets (e.g. TSI)
    discret_->SetState(0,"residual displacement", zeros_);
    discret_->SetState(0,"displacement", (*dis_)(0));
    if (damping_ == INPAR::STR::damp_material) discret_->SetState(0,"velocity", (*vel_)(0));

    //set coupling state for volume coupled problems (e.g. for tsi and poro)
    SetCouplingState();

    discret_->Evaluate(p, stiff_, mass_, fint, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Add(*stiff_, false, dampk_, 0.0);
    damp_->Add(*mass_, false, dampm_, 1.0);
    damp_->Complete();
  }

  // in case of C0 pressure field, we need to get rid of
  // pressure equations
  Teuchos::RCP<LINALG::SparseOperator> mass = Teuchos::null;
  if (pressure_ != Teuchos::null)
  {
    mass = Teuchos::rcp(new LINALG::SparseMatrix(*MassMatrix(),Copy));
    mass->ApplyDirichlet(*(pressure_->CondMap()));
  }
  else if (porositysplitter_ != Teuchos::null)
  {
    mass = Teuchos::rcp(new LINALG::SparseMatrix(*MassMatrix(),Copy));
    mass->ApplyDirichlet(*(porositysplitter_->CondMap()));
  }
  else
  {
    mass = mass_;
  }

  // calculate consistent initial accelerations
  // WE MISS:
  //   - surface stress forces
  //   - potential forces
  {
    Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap_, true);
    if (damping_ == INPAR::STR::damp_rayleigh)
    {
      damp_->Multiply(false, (*vel_)[0], *rhs);
    }
    rhs->Update(-1.0, *fint, 1.0, *fext, -1.0);

    // blank RHS on DBC DOFs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), rhs);
    if (pressure_ != Teuchos::null)
      pressure_->InsertCondVector(pressure_->ExtractCondVector(zeros_), rhs);
    if (porositysplitter_ != Teuchos::null)
      porositysplitter_->InsertCondVector(porositysplitter_->ExtractCondVector(zeros_), rhs);
    solver_->Solve(mass->EpetraOperator(), (*acc_)(0), rhs, true, true);
  }

  // We need to reset the stiffness matrix because its graph (topology)
  // is not finished yet in case of constraints and posssibly other side
  // effects (basically managers).
  stiff_->Reset();

  // leave this hell
  return;
}

/*---------------------------------------------------------------*/
/* Apply Dirichlet boundary conditions on provided state vectors */
void STR::TimInt::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc,
  bool recreatemap
)
{
  // In the case of local coordinate systems, we have to rotate forward ...
  // --------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (dis != Teuchos::null)
        locsysman_->RotateGlobalToLocal(dis,true);
    if (vel != Teuchos::null)
        locsysman_->RotateGlobalToLocal(vel);
    if (acc != Teuchos::null)
        locsysman_->RotateGlobalToLocal(acc);
  }

  // Apply DBCs
  // --------------------------------------------------------------------------------
  // needed parameters
  ParameterList p;
  p.set("total time", time);  // target time

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_->ClearState();
  if (recreatemap)
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                                Teuchos::null, dbcmaps_);
  }
  else
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                               Teuchos::null, Teuchos::null);
  }
  discret_->ClearState();

  // In the case of local coordinate systems, we have to rotate back into global Cartesian frame
  // --------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (dis != Teuchos::null)
        locsysman_->RotateLocalToGlobal(dis,true);
    if (vel != Teuchos::null)
        locsysman_->RotateLocalToGlobal(vel);
    if (acc != Teuchos::null)
        locsysman_->RotateLocalToGlobal(acc);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Update time and step counter */
void STR::TimInt::UpdateStepTime()
{
  // update time and step
  time_->UpdateSteps(timen_);  // t_{n} := t_{n+1}, etc
  step_ = stepn_;  // n := n+1
  //
  timen_ += (*dt_)[0];
  stepn_ += 1;

  // new deal
  return;
}

/*----------------------------------------------------------------------*/
/* Update contact and meshtying */
void STR::TimInt::UpdateStepContactMeshtying()
{
   if (HaveContactMeshtying())
     cmtman_->GetStrategy().Update(stepn_,disn_);

   // ciao
   return;
}

/*----------------------------------------------------------------------*/
/* Update beam contact */
void STR::TimInt::UpdateStepBeamContact()
{
   if (HaveBeamContact())
     beamcman_->Update(*disn_,stepn_,99);

   // ciao
   return;
}

/*----------------------------------------------------------------------*/
/* Velocity update method (VUM) for contact */
void STR::TimInt::UpdateStepContactVUM()
{
  if (HaveContactMeshtying())
  {
    bool do_vum = DRT::INPUT::IntegralValue<int>(cmtman_->GetStrategy().Params(),"VELOCITY_UPDATE");

    //********************************************************************
    // VELOCITY UPDATE METHOD
    //********************************************************************
    if (do_vum)
    {
      // not yet implemented
      dserror("ERROR: Velocity update method not yet implemented");
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Reset configuration after time step */
void STR::TimInt::ResetStep()
{
  // reset state vectors
  disn_->Update(1.0, (*dis_)[0], 0.0);
  if (dismatn_ != Teuchos::null)
    dismatn_->Update(1.0, (*dism_)[0], 0.0);
  veln_->Update(1.0, (*vel_)[0], 0.0);
  accn_->Update(1.0, (*acc_)[0], 0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_struct_reset_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // I am gone
  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart values */
void STR::TimInt::ReadRestart
(
  const int step
)
{
  IO::DiscretizationReader reader(discret_, step);
  if (step != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");

  step_ = step;
  stepn_ = step_ + 1;
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  timen_ = (*time_)[0] + (*dt_)[0];

  ReadRestartState();
  ReadRestartConstraint();
  ReadRestartContactMeshtying();
  ReadRestartStatMech();
  ReadRestartSurfstress();
  ReadRestartMultiScale();

  ReadRestartForce();

  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_->DofRowMap();
}

/*----------------------------------------------------------------------*/
/* Set restart values passed down from above */
void STR::TimInt::SetRestart
(
  int step,
  double time,
  Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector> veln,
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<std::vector<char> > elementdata
)
{
  step_ = step;
  stepn_ = step_ + 1;
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, time));
  timen_ = (*time_)[0] + (*dt_)[0];

  SetRestartState(disn,veln,accn,elementdata);

  // set restart is only for simple structure problems
  // hence we put some security measures in place
  // surface stress

  if (surfstressman_->HaveSurfStress()) dserror("Set restart not implemented for surface stress");

  // constraints
  if (conman_->HaveConstraint()) dserror("Set restart not implemented for constraints");

  // contact / meshtying
  if (HaveContactMeshtying()) dserror("Set restart not implemented for contact / meshtying");

  // statistical mechanics
  if (HaveStatMech()) dserror("Set restart not implemented forstatistical mechanics");

  //biofilm growth
  if (strgrdisp_!=Teuchos::null) dserror("Set restart not implemented for biofilm growth");

  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_->DofRowMap();
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void STR::TimInt::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  reader.ReadVector(disn_, "displacement");
  dis_->UpdateSteps(*disn_);

  if( (dismatn_!=Teuchos::null) )
  {
    reader.ReadVector(dismatn_, "material_displacement");
    dism_->UpdateSteps(*dismatn_);
  }

  reader.ReadVector(veln_, "velocity");
  vel_->UpdateSteps(*veln_);
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);
  reader.ReadMesh(step_);

  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void STR::TimInt::SetRestartState
(
    Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::RCP<Epetra_Vector> accn,
    Teuchos::RCP<std::vector<char> > elementdata
    )
{
  dis_->UpdateSteps(*disn);
  vel_->UpdateSteps(*veln);
  acc_->UpdateSteps(*accn);

  // the following is copied from Readmesh()
  // before we unpack nodes/elements we store a copy of the nodal row/col map
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*discret_->NodeRowMap()));
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(*discret_->NodeColMap()));

  // unpack nodes and elements 
  discret_->UnPackMyElements(elementdata);
  int err =  discret_->FillComplete(true,true,true);
  if (err) dserror("FillComplete() returned err=%d",err);
  return;
}
/*----------------------------------------------------------------------*/
/* Read and set restart values for constraints */
void STR::TimInt::ReadRestartConstraint()
{
  if (conman_->HaveConstraint())
  {
    IO::DiscretizationReader reader(discret_, step_);
    double uzawatemp = reader.ReadDouble("uzawaparameter");
    consolv_->SetUzawaParameter(uzawatemp);
    Teuchos::RCP<Epetra_Map> constrmap=conman_->GetConstraintMap();
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*constrmap, true);
    reader.ReadVector(tempvec, "lagrmultiplier");
    conman_->SetLagrMultVector(tempvec);
    reader.ReadVector(tempvec, "refconval");
    conman_->SetRefBaseValues(tempvec, (*time_)[0]);

  }
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for contact / meshtying */
void STR::TimInt::ReadRestartContactMeshtying()
{
  if (HaveContactMeshtying())
  {
    //**********************************************************************
    // NOTE: There is an important difference here between contact and
    // meshtying simulations. In both cases, the current coupling operators
    // have to be re-computed for restart, but in the meshtying case this
    // means evaluating DM in the reference configuration!
    // Thus, both dis_ (current displacement state) and zero_ are handed
    // in and contact / meshtying managers choose the correct state.
    //**********************************************************************
    IO::DiscretizationReader reader(discret_,step_);
    cmtman_->ReadRestart(reader,(*dis_)(0),zeros_);
  }
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for statmech */
void STR::TimInt::ReadRestartStatMech()
{
  if (HaveStatMech())
  {
    IO::DiscretizationReader reader(discret_,step_);
    statmechman_->ReadRestart(reader, (*dt_)[0]);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for constraints */
void STR::TimInt::ReadRestartSurfstress()
{
  if (surfstressman_ -> HaveSurfStress())
    surfstressman_->ReadRestart(step_, DRT::Problem::Instance()->InputControlFile()->FileName());
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for multi-scale */
void STR::TimInt::ReadRestartMultiScale()
{
  Teuchos::RCP<MAT::PAR::Bundle> materials = DRT::Problem::Instance()->Materials();
  for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials->Map()->begin();
       i!=materials->Map()->end();
       ++i)
  {
    Teuchos::RCP<MAT::PAR::Material> mat = i->second;
    if (mat->Type() == INPAR::MAT::m_struct_multiscale)
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action", "multi_readrestart");
      discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                         Teuchos::null, Teuchos::null, Teuchos::null);
      discret_->ClearState();
      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/* Calculate all output quantities that depend on a potential material history */
void STR::TimInt::PrepareOutput()
{
  DetermineStressStrain();
  DetermineEnergy();
  if (havemicromat_) PrepareOutputMicro();
}

/*----------------------------------------------------------------------*/
/* output to file
 * originally by mwgee 03/07 */
void STR::TimInt::OutputStep(bool forced_writerestart)
{
  // special treatment is necessary when restart is forced
  if(forced_writerestart)
  {
    // reset possible history data on element level
    ResetStep();
    // restart has already been written
    if(writerestartevery_ and (step_%writerestartevery_ == 0))
      return;
    // if state already exists, add restart information
    if(writeresultsevery_ and (step_%writeresultsevery_ == 0) and step_!=DRT::Problem::Instance()->Restart())
    {
      AddRestartToOutputState();
      return;
    }
  }

  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if ( (writerestartevery_ and (step_%writerestartevery_ == 0)) or forced_writerestart )
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if ( writestate_
       and writeresultsevery_ and (step_%writeresultsevery_ == 0)
       and (not datawritten) )
  {
    OutputState(datawritten);
  }

  // output stress & strain
  if ( writeresultsevery_
       and ( (writestress_ != INPAR::STR::stress_none)
             or (writecouplstress_ != INPAR::STR::stress_none)
             or (writestrain_ != INPAR::STR::strain_none)
             or (writeplstrain_ != INPAR::STR::strain_none) )
       and (step_%writeresultsevery_ == 0) )
  {
    OutputStressStrain(datawritten);
  }

  // output energy
  if ( writeenergyevery_ and (step_%writeenergyevery_ == 0) )
  {
    OutputEnergy();
  }

  // output active set, energies and momentum for contact
  OutputContact();

  // print error norms
  OutputErrorNorms();

  // output of nodal positions in current configuration
  OutputNodalPositions();

  // write output on micro-scale (multi-scale analysis)
  if (havemicromat_) OutputMicro();

  // write patient specific output
  if (writeresultsevery_ and (step_%writeresultsevery_ == 0))
  {
    OutputPatspec();
  }

  // what's next?
  return;
}
/*----------------------------------------------------------------------*/
/* We need the restart data to perform on "restarts" on the fly for parameter
 * continuation
 */
void STR::TimInt::GetRestartData
(
  Teuchos::RCP<int> step,
  Teuchos::RCP<double> time,
  Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector> veln,
  Teuchos::RCP<Epetra_Vector> accn,
  Teuchos::RCP<std::vector<char> > elementdata
)
{
  // at some point we have to create a copy
  *step = step_ ;
  *time = (*time_)[0];
  *disn = *disn_;
  *veln = *veln_;
  *accn = *accn_;
  *elementdata = *(discret_->PackMyElements());

  // get restart data is only for simple structure problems
  // hence

  // surface stress
  if (surfstressman_->HaveSurfStress()) dserror("Get restart data not implemented for surface stress");

  // constraints
  if (conman_->HaveConstraint()) dserror("Get restart data not implemented for constraints");

  // contact / meshtying
  if (HaveContactMeshtying()) dserror("Get restart data not implemented for contact / meshtying");

  // statistical mechanics
  if (HaveStatMech()) dserror("Get restart data not implemented for statistical mechanics");

  //biofilm growth
  if (strgrdisp_!=Teuchos::null) dserror("Get restart data not implemented for biofilm growth");

  return;
}
/*----------------------------------------------------------------------*/
/* write restart
 * originally by mwgee 03/07 */
void STR::TimInt::OutputRestart
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;
  // for multilevel monte carlo we do not need to write mesh in every run
  if (writereducedrestart_==true)
  {
    // write restart output, please
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    output_->WriteElementData(firstoutputofrun_);
  }
  else
  {
    // write restart output, please
    output_->WriteMesh(step_, (*time_)[0]);
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    if( dism_!=Teuchos::null )
      output_->WriteVector("material_displacement", (*dism_)(0));
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
    if(!HaveStatMech())
      output_->WriteElementData(firstoutputofrun_);
    WriteRestartForce(output_);
  }
  // owner of elements is just written once because it does not change during simulation (so far)
  firstoutputofrun_ = false;

  // surface stress
  if (surfstressman_->HaveSurfStress())
  {
    surfstressman_->WriteRestart(step_, (*time_)[0]);
  }

  // constraints
  if (conman_->HaveConstraint())
  {
    output_->WriteDouble("uzawaparameter",
                          consolv_->GetUzawaParameter());
    output_->WriteVector("lagrmultiplier",
                          conman_->GetLagrMultVector());
    output_->WriteVector("refconval",
                          conman_->GetRefBaseValues());
  }

  // contact / meshtying
  if (HaveContactMeshtying())
  {
      cmtman_->WriteRestart(*output_);

      cmtman_->PostprocessTractions(*output_);
  }

  // statistical mechanics
  if (HaveStatMech())
  {
    statmechman_->WriteRestart(output_, (*dt_)[0]);
  }

  //biofilm growth
  if (strgrdisp_!=Teuchos::null)
  {
    output_->WriteVector("str_growth_displ", strgrdisp_);
  }

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0))
  {
    IO::cout << "====== Restart written in step " << step_ << IO::endl;
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart written in step %d\n", step_);
    fflush(errfile_);
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations
 * originally by mwgee 03/07 */
void STR::TimInt::OutputState
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // write now
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", (*dis_)(0));

  if( (dismatn_!=Teuchos::null))
    output_->WriteVector("material_displacement", (*dism_)(0));

  // for visualization of vel and acc do not forget to comment in corresponding lines in StructureEnsightWriter
  if(writevelacc_)
  {
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
  }

  //biofilm growth
  if (strgrdisp_!=Teuchos::null)
  {
    output_->WriteVector("str_growth_displ", strgrdisp_);
  }

  // owner of elements is just written once because it does not change during simulation (so far)
  output_->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  if (surfstressman_->HaveSurfStress() && writesurfactant_)
    surfstressman_->WriteResults(step_, (*time_)[0]);

  // contact / meshtying
  if (HaveContactMeshtying())
    cmtman_->PostprocessTractions(*output_);

  if(porositysplitter_!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> porosity = porositysplitter_->ExtractCondVector((*dis_)(0));
    output_->WriteVector("porosity_p1", porosity);
  }
  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* add restart information to OutputState */
void STR::TimInt::AddRestartToOutputState()
{
  // add velocity and acceleration if necessary
  if(!writevelacc_)
  {
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
  }

  WriteRestartForce(output_);

  // constraints
  if (conman_->HaveConstraint())
  {
    output_->WriteDouble("uzawaparameter",
                          consolv_->GetUzawaParameter());
    output_->WriteVector("lagrmultiplier",
                          conman_->GetLagrMultVector());
    output_->WriteVector("refconval",
                          conman_->GetRefBaseValues());
  }

  // contact/meshtying information
  if (HaveContactMeshtying())
  {
    cmtman_->WriteRestart(*output_,true);
  }

  // TODO: add missing restart data for surface stress, contact/meshtying and StatMech here


  // finally add the missing mesh information, order is important here
  output_->WriteMesh(step_, (*time_)[0]);

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0))
  {
    IO::cout << "====== Restart written in step " << step_ << IO::endl;
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart written in step %d\n", step_);
    fflush(errfile_);
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* Calculation of stresses and strains */
void STR::TimInt::DetermineStressStrain()
{
  if ( writeresultsevery_
       and ( (writestress_ != INPAR::STR::stress_none)
             or (writecouplstress_ != INPAR::STR::stress_none)
             or (writestrain_ != INPAR::STR::strain_none)
             or (writeplstrain_ != INPAR::STR::strain_none) )
       and (stepn_%writeresultsevery_ == 0) )
  {
    //-------------------------------
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", (*dt_)[0]);
    p.set<int>("young_temp", young_temp_);

    stressdata_ = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stressdata_);
    p.set<int>("iostress", writestress_);

    // write stress data that arise from the coupling with another field, e.g.
    // in TSI: couplstress corresponds to thermal stresses
    couplstressdata_ = Teuchos::rcp(new std::vector<char>());
    p.set("couplstress", couplstressdata_);
    p.set<int>("iocouplstress", writecouplstress_);

    straindata_ = Teuchos::rcp(new std::vector<char>());
    p.set("strain", straindata_);
    p.set<int>("iostrain", writestrain_);

    // plastic strain
    plstraindata_ = Teuchos::rcp(new std::vector<char>());
    p.set("plstrain", plstraindata_);
    p.set<int>("ioplstrain", writeplstrain_);

    // set plasticity data
    if (HaveSemiSmoothPlasticity()) plastman_->SetPlasticParams(p);

    // set vector values needed by elements
    discret_->ClearState();
    // extended SetState(0,...) in case of multiple dofsets (e.g. TSI)
    discret_->SetState(0,"residual displacement", zeros_);
    discret_->SetState(0,"displacement", disn_);

    //set coupling state for volume coupled problems (e.g. for tsi and poro)
    SetCouplingState();

    if( (dismatn_!=Teuchos::null) )
      discret_->SetState(0,"material_displacement",dismatn_);

    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();

    // set plasticity data
    if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(p);
  }
}

/*----------------------------------------------------------------------*/
/* Calculation of internal, external and kinetic energy */
void STR::TimInt::DetermineEnergy()
{
  if ( writeenergyevery_ and (stepn_%writeenergyevery_ == 0) )
  {
    // internal/strain energy
    intergy_ = 0.0;  // total internal energy
    {
      ParameterList p;
      // other parameters needed by the elements
      p.set("action", "calc_struct_energy");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("displacement", disn_);
      // get energies
      Teuchos::RCP<Epetra_SerialDenseVector> energies
        = Teuchos::rcp(new Epetra_SerialDenseVector(1));
      discret_->EvaluateScalars(p, energies);
      discret_->ClearState();
      intergy_ = (*energies)(0);
    }

    // global calculation of kinetic energy
    kinergy_ = 0.0;  // total kinetic energy
    {
      Teuchos::RCP<Epetra_Vector> linmom
        = LINALG::CreateVector(*dofrowmap_, true);
      mass_->Multiply(false, *veln_, *linmom);
      linmom->Dot(*veln_, &kinergy_);
      kinergy_ *= 0.5;
    }

    // external energy
    extergy_ = 0.0;  // total external energy
    {
      // WARNING: This will only work with dead loads and implicit time
      // integration (otherwise there is no fextn_)!!!
      Teuchos::RCP<Epetra_Vector> fext = FextNew();
      fext->Dot(*disn_, &extergy_);
    }
  }
}

/*----------------------------------------------------------------------*/
/* stress calculation and output */
void STR::TimInt::OutputStressStrain
(
  bool& datawritten
)
{
  // Make new step
  if (not datawritten)
  {
    output_->NewStep(step_, (*time_)[0]);
  }
  datawritten = true;

  // write stress
  if (writestress_ != INPAR::STR::stress_none)
  {
    std::string stresstext = "";
    if (writestress_ == INPAR::STR::stress_cauchy)
    {
      stresstext = "gauss_cauchy_stresses_xyz";
    }
    else if (writestress_ == INPAR::STR::stress_2pk)
    {
      stresstext = "gauss_2PK_stresses_xyz";
    }
    else
    {
      dserror("requested stress type not supported");
    }
    output_->WriteVector(stresstext, *stressdata_,
                         *(discret_->ElementRowMap()));
    // we don't need this anymore
    stressdata_ = Teuchos::null;
  }

  // write coupling stress
  if (writecouplstress_ != INPAR::STR::stress_none)
  {
    std::string couplstresstext = "";
    if (writecouplstress_ == INPAR::STR::stress_cauchy)
    {
      couplstresstext = "gauss_cauchy_coupling_stresses_xyz";
    }
    else if (writecouplstress_ == INPAR::STR::stress_2pk)
    {
      couplstresstext = "gauss_2PK_coupling_stresses_xyz";
    }
    else
    {
      dserror("requested stress type not supported");
    }
    output_->WriteVector(couplstresstext, *couplstressdata_,
                         *(discret_->ElementRowMap()));
    // we don't need this anymore
    couplstressdata_ = Teuchos::null;
  }

  // write strain
  if (writestrain_ != INPAR::STR::strain_none)
  {
    std::string straintext = "";
    if (writestrain_ == INPAR::STR::strain_ea)
    {
      straintext = "gauss_EA_strains_xyz";
    }
    else if (writestrain_ == INPAR::STR::strain_gl)
    {
      straintext = "gauss_GL_strains_xyz";
    }
    else if (writestrain_ == INPAR::STR::strain_log)
    {
      straintext = "gauss_LOG_strains_xyz";
    }
    else
    {
      dserror("requested strain type not supported");
    }
    output_->WriteVector(straintext, *straindata_,
                         *(discret_->ElementRowMap()));
    // we don't need this anymore
    straindata_ = Teuchos::null;
  }

  // write plastic strain
  if (writeplstrain_ != INPAR::STR::strain_none)
  {
    std::string plstraintext = "";
    if (writeplstrain_ == INPAR::STR::strain_ea)
    {
      plstraintext = "gauss_pl_EA_strains_xyz";
    }
    else if (writeplstrain_ == INPAR::STR::strain_gl)
    {
      plstraintext = "gauss_pl_GL_strains_xyz";
    }
    else
    {
      dserror("requested plastic strain type not supported");
    }
    output_->WriteVector(plstraintext, *plstraindata_,
                         *(discret_->ElementRowMap()));
    // we don't need this anymore
    plstraindata_ = Teuchos::null;
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
/* output system energies */
void STR::TimInt::OutputEnergy()
{
  // total energy
  double totergy = kinergy_ + intergy_ - extergy_;

  // the output
  if (myrank_ == 0)
  {
    *energyfile_ << " " << std::setw(9) << step_
                 << std::scientific  << std::setprecision(16)
                 << " " << (*time_)[0]
                 << " " << totergy
                 << " " << kinergy_
                 << " " << intergy_
                 << " " << extergy_
                 << std::endl;
  }

  // in God we trust
  return;
}

/*----------------------------------------------------------------------*/
/* output active set, energies and momentum for contact */
void STR::TimInt::OutputContact()
{
  // only for contact / meshtying simulations
  if (HaveContactMeshtying())
  {
    // THIS IS FOR DEBUGGING ONLY!!!
    // print contact forces with respect to reference configuration
  #ifdef CONTACTFORCEREFCONFIG
    cmtman_->GetStrategy().ForceRefConfig();
  #endif

    // print active set
    cmtman_->GetStrategy().PrintActiveSet();

    // check chosen output option
    INPAR::CONTACT::EmOutputType emtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::EmOutputType>(cmtman_->GetStrategy().Params(),"EMOUTPUT");

    // get out of here if no enrgy momentum output wanted
    if (emtype==INPAR::CONTACT::output_none) return;

    // get some parameters from parameter list
    double timen = (*time_)[0];
    double dt    = (*dt_)[0];
    int dim      = cmtman_->GetStrategy().Dim();

    // global linear momentum (M*v)
    RCP<Epetra_Vector> mv = LINALG::CreateVector(*(discret_->DofRowMap()), true);
    mass_->Multiply(false, (*vel_)[0], *mv);

    // linear / angular momentum
    std::vector<double> sumlinmom(3);
    std::vector<double> sumangmom(3);
    std::vector<double> angmom(3);
    std::vector<double> linmom(3);

    // vectors of nodal properties
    std::vector<double> nodelinmom(3);
    std::vector<double> nodeangmom(3);
    std::vector<double> position(3);

    // loop over all nodes belonging to the respective processor
    for (int k=0; k<(discret_->NodeRowMap())->NumMyElements();++k)
    {
      // get current node
      int gid = (discret_->NodeRowMap())->GID(k);
      DRT::Node* mynode = discret_->gNode(gid);
      std::vector<int> globaldofs = discret_->Dof(mynode);

      // loop over all DOFs comprised by this node
      for (int i=0;i<dim;i++)
      {
        nodelinmom[i] = (*mv)[mv->Map().LID(globaldofs[i])];
        sumlinmom[i] += nodelinmom[i];
        position[i]   = (mynode->X())[i] + ((*dis_)[0])[mv->Map().LID(globaldofs[i])];
      }

      // calculate vector product position x linmom
      nodeangmom[0] = position[1]*nodelinmom[2] - position[2]*nodelinmom[1];
      nodeangmom[1] = position[2]*nodelinmom[0] - position[0]*nodelinmom[2];
      nodeangmom[2] = position[0]*nodelinmom[1] - position[1]*nodelinmom[0];

      // loop over all DOFs comprised by this node
      for (int i=0; i<3; ++i) sumangmom[i] += nodeangmom[i];
    }

    // global quantities (sum over all processors)
    for (int i=0;i<3;++i)
    {
      cmtman_->Comm().SumAll(&sumangmom[i],&angmom[i],1);
      cmtman_->Comm().SumAll(&sumlinmom[i],&linmom[i],1);
    }

    //--------------------------Calculation of total kinetic energy
    double kinen = 0.0;
    mv->Dot((*vel_)[0],&kinen);
    kinen *= 0.5;

    //-------------------------Calculation of total internal energy
    // BE CAREFUL HERE!!!
    // Currently this is only valid for materials without(!) history
    // data, because we have already updated everything and thus
    // any potential material history has already been overwritten.
    // When we are interested in internal energy output for contact
    // or meshtying simulations with(!) material history, we should
    // move the following block before the first UpdateStep() method.
    double inten = 0.0;
    ParameterList p;
    p.set("action", "calc_struct_energy");
    discret_->ClearState();
    discret_->SetState("displacement", (*dis_)(0));
    RCP<Epetra_SerialDenseVector> energies = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    energies->Scale(0.0);
    discret_->EvaluateScalars(p, energies);
    discret_->ClearState();
    inten = (*energies)(0);

    //-------------------------Calculation of total external energy
    double exten = 0.0;
    // WARNING: This will only work with dead loads!!!
    //fext_->Dot(*dis_, &exten);

    //----------------------------------------Print results to file
    if (emtype == INPAR::CONTACT::output_file ||
        emtype == INPAR::CONTACT::output_both)
    {
      // processor 0 does all the work
      if (!myrank_)
      {
        // path and filename
        std::ostringstream filename;
        const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
        filename << filebase << ".energymomentum";

        // open file
        FILE* MyFile = NULL;
        if (timen < 2*dt)
        {
          MyFile = fopen(filename.str().c_str(),"wt");

          // initialize file pointer for writing contact interface forces/moments
          FILE* MyConForce = NULL;
          std::ostringstream filenameif;
          filenameif << filebase << ".energymomentum";
          MyConForce = fopen(filenameif.str().c_str(), "wt");
          if (MyConForce!=NULL) fclose(MyConForce);
          else dserror("ERROR: File for writing contact interface forces/moments could not be generated.");
        }
        else
          MyFile = fopen(filename.str().c_str(),"at+");

        // add current values to file
        if (MyFile!=NULL)
        {
         std::stringstream filec;
         fprintf(MyFile, "% e\t", timen);
         for (int i=0; i<3; i++) fprintf(MyFile, "% e\t", linmom[i]);
         for (int i=0; i<3; i++) fprintf(MyFile, "% e\t", angmom[i]);
         fprintf(MyFile, "% e\t% e\t% e\t% e\n",kinen,inten,exten,kinen+inten-exten);
         fclose(MyFile);
        }
        else
          dserror("ERROR: File for writing momentum and energy data could not be opened.");
      }
    }

    //-------------------------------Print energy results to screen
    if (emtype == INPAR::CONTACT::output_screen ||
        emtype == INPAR::CONTACT::output_both)
    {
      // processor 0 does all the work
      if (!myrank_)
      {
        printf("******************************");
        printf("\nMECHANICAL ENERGIES:");
        printf("\nE_kinetic \t %e",kinen);
        printf("\nE_internal \t %e",inten);
        printf("\nE_external \t %e",exten);
        printf("\n------------------------------");
        printf("\nE_total \t %e",kinen+inten-exten);
        printf("\n******************************");

        printf("\n\n********************************************");
        printf("\nLINEAR / ANGULAR MOMENTUM:");
        printf("\nL_x  % e \t H_x  % e",linmom[0],angmom[0]);
        printf("\nL_y  % e \t H_y  % e",linmom[1],angmom[1]);
        printf("\nL_z  % e \t H_z  % e",linmom[2],angmom[2]);
        printf("\n********************************************\n\n");
        fflush(stdout);
      }
    }

    //-------------------------- Compute and output interface forces
    cmtman_->GetStrategy().InterfaceForces(true);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* output solution error norms */
void STR::TimInt::OutputErrorNorms()
{
  // get out of here if no output wanted
  const Teuchos::ParameterList& listcmt = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::CONTACT::ErrorNorms entype = DRT::INPUT::IntegralValue<INPAR::CONTACT::ErrorNorms>(listcmt,"ERROR_NORMS");
  if (entype==INPAR::CONTACT::errornorms_none) return;

  // initialize variables
  RCP<Epetra_SerialDenseVector> norms = Teuchos::rcp(new Epetra_SerialDenseVector(3));
  norms->Scale(0.0);

  // call discretization to evaluate error norms
  ParameterList p;
  p.set("action", "calc_struct_errornorms");
  discret_->ClearState();
  discret_->SetState("displacement",(*dis_)(0));
  discret_->EvaluateScalars(p, norms);
  discret_->ClearState();

  // proc 0 writes output to screen
  if (!myrank_)
  {
    printf("**********************************");
    printf("\nSOLUTION ERROR NORMS:");
    printf("\nL_2 norm:     %.10e",sqrt((*norms)(0)));
    printf("\nH_1 norm:     %.10e",sqrt((*norms)(1)));
    printf("\nEnergy norm:  %.10e",sqrt((*norms)(2)));
    printf("\n**********************************\n\n");
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* output on micro-scale */
void STR::TimInt::OutputMicro()
{
  for (int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele = discret_->lRowElement(i);
    RCP<MAT::Material> mat = actele->Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Output();
    }
  }
}

/*----------------------------------------------------------------------*/
/* calculate stresses and strains on micro-scale */
void STR::TimInt::PrepareOutputMicro()
{
  for (int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele = discret_->lRowElement(i);

    RCP<MAT::Material> mat = actele->Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->PrepareOutput();
    }
  }
}

/*----------------------------------------------------------------------*/
/* output nodal positions */
void STR::TimInt::OutputNodalPositions()
{
#ifdef PRINTSTRUCTDEFORMEDNODECOORDS

  /////////////////////////////////////////////////////////////////
  // from here I want to output my ale displacements - devaal 14.12.2010
  /////////////////////////////////////////////////////////////////

  if (discret_->Comm().NumProc() != 1)
    dserror("The flag PRINTSTRUCTDEFORMEDNODECOORDS is on and only works with 1 processor");

  std::cout << "STRUCT DISCRETIZATION IN THE DEFORMED CONFIGURATIONS" << std::endl;
  // does discret_ exist here?
  //std::cout << "discret_->NodeRowMap()" << discret_->NodeRowMap() << std::endl;

  //RCP<Epetra_Vector> mynoderowmap = Teuchos::rcp(new Epetra_Vector(discret_->NodeRowMap()));
  //RCP<Epetra_Vector> noderowmap_ = Teuchos::rcp(new Epetra_Vector(discret_->NodeRowMap()));
  //dofrowmap_  = Teuchos::rcp(new discret_->DofRowMap());
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  const Epetra_Map* dofrowmap = discret_->DofRowMap();


  for (int lid=0; lid<noderowmap->NumGlobalPoints(); lid++)
  {
    int gid;
    // get global id of a node
    gid = noderowmap->GID(lid);
    // get the node
    DRT::Node * node = discret_->gNode(gid);
    //std::cout<<"mynode"<<*node<<std::endl;

    // get the coordinates of the node
    const double * X = node->X();
    // get degrees of freedom of a node

    std::vector<int> gdofs = discret_->Dof(node);
    //std::cout << "for node:" << *node << std::endl;
    //std::cout << "this is my gdof vector" << gdofs[0] << " " << gdofs[1] << " " << gdofs[2] << std::endl;

    // get displacements of a node
    std::vector<double> mydisp (3,0.0);
    for (int ldof = 0; ldof<3; ldof ++)
    {
      int displid = dofrowmap->LID(gdofs[ldof]);

      //std::cout << "displacement local id - in the rowmap" << displid << std::endl;
      mydisp[ldof] = ((*dis_)[0])[displid];
      //std::cout << "at node" << gid << "mydisplacement in each direction" << mydisp[ldof] << std::endl;
      //make zero if it is too small
      if (abs(mydisp[ldof]) < 0.00001)
      {
        mydisp[ldof] = 0.0;
      }
    }

    // Export disp, X
    double newX = mydisp[0]+X[0];
    double newY = mydisp[1]+X[1];
    double newZ = mydisp[2]+X[2];
    //std::cout << "NODE " << gid << "  COORD  " << newX << " " << newY << " " << newZ << std::endl;
    std::cout << gid << " " << newX << " " << newY << " " << newZ << std::endl;
  }

#endif //PRINTSTRUCTDEFORMEDNODECOORDS

  return;
}

/*----------------------------------------------------------------------*/
/* output patient specific stuff */
void STR::TimInt::OutputPatspec()
{
  // do the output for the patient specific conditions (if they exist)
  const Teuchos::ParameterList& patspec  = DRT::Problem::Instance()->PatSpecParams();
  if (DRT::INPUT::IntegralValue<int>(patspec,"PATSPEC"))
  {
    PATSPEC::PatspecOutput(output_,discret_,pslist_);
  }
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void STR::TimInt::ApplyForceExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< old displacement state
  const Teuchos::RCP<Epetra_Vector> disn,  //!< new displacement state
  const Teuchos::RCP<Epetra_Vector> vel,  //!< velocity state
  Teuchos::RCP<Epetra_Vector>& fext,  //!< external force
  Teuchos::RCP<LINALG::SparseOperator>& fextlin //!<linearization of external force
)
{
  ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"displacement", dis);

  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);
  // get load vector
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  bool loadlin = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "LOADLIN");
  if (!loadlin)
    discret_->EvaluateNeumann(p, *fext);
  else
  {
    discret_->SetState(0,"displacement new", disn);
    discret_->EvaluateNeumann(p, fext, fextlin);
  }

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at state */
void STR::TimInt::ApplyForceStiffInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // residual displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  Teuchos::RCP<Epetra_Vector> fint,  // internal force
  Teuchos::RCP<LINALG::SparseOperator> stiff,  // stiffness matrix
  Teuchos::RCP<LINALG::SparseOperator> damp   // material damping matrix
)
{
  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  const std::string action = "calc_struct_nlnstiff";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  p.set("damping", damping_);
  p.set<int>("young_temp", young_temp_);
  if (pressure_ != Teuchos::null) p.set("volume", 0.0);

  // set plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->SetPlasticParams(p);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"residual displacement", disi);
  discret_->SetState(0,"displacement", dis);
  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector

  //set coupling state for volume coupled problems (e.g. for tsi and poro)
  SetCouplingState();

  // Set material displacement state for ale-wear formulation
  if( (dismatn_!=Teuchos::null) )
    discret_->SetState(0,"material_displacement",dismatn_);

  discret_->Evaluate(p, stiff, damp, fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // get plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(p);

#if 0
  if (pressure_ != Teuchos::null)
    std::cout << "Total volume=" << std::scientific << p.get<double>("volume") << std::endl;
#endif

  // *********** time measurement ***********
  dtele_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate inertial force and its linearization */
void STR::TimInt::ApplyForceStiffInternalAndInertial
(
  const double time,  //!< evaluation time
  const double dt,  //!< step size
  const double timintfac_dis,  //!< additional parameter from time integration 1
  const double timintfac_vel,  //!< additional parameter from time integration 2
  const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  const Teuchos::RCP<Epetra_Vector> acc,  // acceleration state
  Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
  Teuchos::RCP<Epetra_Vector> finert,  //!< inertial force
  Teuchos::RCP<LINALG::SparseOperator> stiff,  //!< stiffness matrix
  Teuchos::RCP<LINALG::SparseOperator> mass  //!< mass matrix
)
{
    //dis->Print(std::cout);
    //vel->Print(std::cout);
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // action for elements
    const std::string action = "calc_struct_nlnstiffmass";
    p.set("action", action);
    // other parameters that might be needed by the elements
    p.set("total time", time);
    p.set("delta time", dt);

    p.set("timintfac_dis", timintfac_dis);
    p.set("timintfac_vel", timintfac_vel);

    // set plasticity data
    if (HaveSemiSmoothPlasticity()) plastman_->SetPlasticParams(p);

    discret_->ClearState();
    discret_->SetState(0,"residual displacement", disi);
    discret_->SetState(0,"displacement", dis);
    discret_->SetState(0,"velocity", vel);
    discret_->SetState(0,"acceleration", acc);

    discret_->Evaluate(p, stiff, mass, fint, finert, Teuchos::null);
    discret_->ClearState();

    // get plasticity data
    if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(p);

    mass->Complete();

  return;
};

/*----------------------------------------------------------------------*/
/* check wether we have nonlinear inertia forces or not */
bool STR::TimInt::HaveNonlinearMass()
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  bool masslin = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "MASSLIN");

  return masslin;
}

/*----------------------------------------------------------------------*/
/* check wether the initial conditions are fulfilled */
void STR::TimInt::NonlinearMassSanityCheck( Teuchos::RCP<Epetra_Vector> fext,
                                            Teuchos::RCP<Epetra_Vector> dis,
                                            const Teuchos::RCP<Epetra_Vector> vel,
                                            const Teuchos::RCP<Epetra_Vector> acc)
{
  double fextnorm;
  fext->Norm2(&fextnorm);

  double dispnorm;
  dis->Norm2(&dispnorm);

  double velnorm;
  vel->Norm2(&velnorm);

  double accnorm;
  acc->Norm2(&accnorm);

  if (fextnorm > 1.0e-14)
    dserror("Initial configuration does not fulfill equilibrium, check your initial external forces, velocities and accelerations!!!");

  if ((dispnorm > 1.0e-14) or (velnorm > 1.0e-14) or (accnorm > 1.0e-14))
    dserror("Nonlinear inertia terms (input flag 'MASSLIN' set to 'yes') are only possible for vanishing initial displacements, velocities and accelerations so far!!!");

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force */
void STR::TimInt::ApplyForceInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // incremental displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  Teuchos::RCP<Epetra_Vector> fint  // internal force
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  std::string action = "calc_struct_internalforce";

  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  p.set<int>("young_temp", young_temp_);

  // set plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->SetPlasticParams(p);

  if (pressure_ != Teuchos::null) p.set("volume", 0.0);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement", disi);  // these are incremental
  discret_->SetState("displacement", dis);

  //set coupling state for volume coupled problems (e.g. for tsi and poro)
  SetCouplingState();

  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     fint, Teuchos::null, Teuchos::null);

  // get plasticity data
  if (HaveSemiSmoothPlasticity()) plastman_->GetPlasticParams(p);

  discret_->ClearState();

  // where the fun starts
  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
int STR::TimInt::Integrate()
{
  // error checking variables
  int lnonlinsoldiv = 0;
  int nonlinsoldiv  = 0;

  // target time #timen_ and step #stepn_ already set
  // time loop
  // (NOTE: popp 03/2010: we have to add eps to avoid the very
  // awkward effect that the time loop stops one step too early)
  double eps = 1.0e-12;
  while ( (timen_ <= timemax_+eps) and (stepn_ <= stepmax_) and (not nonlinsoldiv) )
  {
    // prepare contact for new time step
    PrepareStepContact();

    // integrate time step
    // after this step we hold disn_, etc
    lnonlinsoldiv = IntegrateStep();
    // since it is possible that the nonlinear solution fails only on some procs
    // we need to communicate the error
    discret_->Comm().Barrier();
    discret_->Comm().MaxAll(&lnonlinsoldiv,&nonlinsoldiv,1);
    discret_->Comm().Barrier();

    // if everything is fine
    if(!nonlinsoldiv)
    {
      // calculate stresses, strains and energies
      // note: this has to be done before the update since otherwise a potential
      // material history is overwritten
      PrepareOutput();

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      UpdateStepState();

      // update time and step
      UpdateStepTime();

      // update everything on the element level
      UpdateStepElement();

      // write output
      OutputStep();

      // print info about finished time step
      PrintStep();
    }
    else // something went wrong update errorcode according to chosen divcont action
    {
      nonlinsoldiv = PerformErrorAction(nonlinsoldiv);
    }
  }
  // stop supporting processors in multi scale simulations
  if (havemicromat_)
  {
    STRUMULTI::stop_np_multiscale();
  }
  // that's it say what went wrong
  return nonlinsoldiv;
}

/*----------------------------------------------------------------------*/
int STR::TimInt::PerformErrorAction(int nonlinsoldiv)
{
  // what to do when nonlinear solver does not converge
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  enum INPAR::STR::DivContAct divcontype = (DRT::INPUT::IntegralValue<INPAR::STR::DivContAct>(sdyn,"DIVERCONT"));
  switch (divcontype)
     {
       case INPAR::STR::divcont_stop:
       {
         // write restart output of last converged step before stopping
         OutputStep(true);

         // we should not get here, dserror for safety
         dserror("Nonlinear solver did not converge! ");
         return 1;
       }
       case INPAR::STR::divcont_continue:
       {
       // we should not get here, dserror for safety
         dserror("Nonlinear solver did not converge! ");
         return 1;
       }
       break;
       case INPAR::STR::divcont_repeat_step:
       {
         IO::cout << "Nonlinear solver failed to converge repeat time step" << IO::endl;
         // do nothing since we didn't update yet
         return 0;
       }
       break;
       case INPAR::STR::divcont_halve_step:
        {
          IO::cout << "Nonlinear solver failed to converge divide timestep in half" << IO::endl;
          // halve the time step size
          (*dt_)[0]=(*dt_)[0]*0.5;
          // update the number of max timesteps
          stepmax_= stepmax_ + (stepmax_-stepn_)*2+1;
          // reset timen_ because it is set in the constructor
          timen_ = (*time_)[0] + (*dt_)[0];;
          return 0;
        }
        break;
       case INPAR::STR::divcont_repeat_simulation:
       {
         if(nonlinsoldiv==1)
           IO::cout << "Nonlinear solver failed to converge and DIVERCONT = repeat_simulation, hence leaving structural time integration " << IO::endl;
         else if (nonlinsoldiv==2)
           IO::cout << "Linear solver failed to converge and DIVERCONT = repeat_simulation, hence leaving structural time integration " << IO::endl;
         return 1; // so that timeloop will be aborted
       }
       break;
       default:
         dserror("Unknown DIVER_CONT case");
         return 1;
         break;
     }
  return 0; // make compiler happy
}
/*----------------------------------------------------------------------*/
/* check whether contact solver should be used */
bool STR::TimInt::UseContactSolver()
{
  // no contact possible -> return false
  if (!HaveContactMeshtying())
    return false;
  // contact possible -> check current status
  else
  {
    // currently not in contact -> return false
    if (!cmtman_->GetStrategy().IsInContact() &&
        !cmtman_->GetStrategy().WasInContact() &&
        !cmtman_->GetStrategy().WasInContactLastTimeStep())
      return false;
    // currently in contact -> return true
    else
      return true;
  }
}

/*----------------------------------------------------------------------*/
/* set volume coupling state from other discretization  vuong 01/12*/
void STR::TimInt::ApplyCouplingState(
  Teuchos::RCP<const Epetra_Vector> state,
  const std::string& name,
  unsigned dofset
  )
{
  //check
  if(state == Teuchos::null)
    dserror("coupling state '%s' is Teuchos::null!",name.c_str());

  if(dofset == 0)
    dserror("dofset number equals zero!");

  //copy state values
  RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(discret_->DofRowMap(dofset)), true);
  LINALG::Export(*state,*tmp);

  //resize the coupling state vector, if it is too short
  if (couplstate_.size()<dofset)
    couplstate_.resize(dofset);

  //save state
  couplstate_[dofset-1][name] = tmp;

  return;
}

/*----------------------------------------------------------------------*/
/* set coupling state for volume coupled problems          */
void STR::TimInt::SetCouplingState()
{
  std::map<std::string,Teuchos::RCP<const Epetra_Vector> >::iterator it;

  //loop over dofsets
  for (unsigned int iter = 0; iter<couplstate_.size();  ++iter)
  {
    //the number of the dof set is iter+1 as couplstate_ does not include the own dof set
    const int dofset = iter+1;

    //set all saved states onto the respective dofset
    for(it=couplstate_[iter].begin(); it!=couplstate_[iter].end(); ++it)
      discret_->SetState(dofset,it->first,it->second);
  }

  return;
}


/*----------------------------------------------------------------------*/
/* Set forces due to interface with fluid,
 * the force is expected external-force-like */
void STR::TimInt::SetForceInterface
(
  Teuchos::RCP<Epetra_MultiVector> iforce  ///< the force on interface
)
{
  fifc_->Update(1.0, *iforce, 0.0);
}

/*----------------------------------------------------------------------*/
/* apply the new material_displacements                      mgit 05/11 */
void STR::TimInt::ApplyDisMat(
  Teuchos::RCP<Epetra_Vector> dismat
  )
{
  // FIXGIT: This is done only for nonzero entries --> try to store the zero-entries to !!!
  // These values are replaced because here, the new absolute material
  // displacement has been evaluated (not the increment)

   for (int k=0;k<dismat->MyLength();++k)
   {
     if ((*dismat)[k] != 0.0)
     {
       (*dismatn_)[k]=(*dismat)[k];
     }
   }
   return;
 }

/*----------------------------------------------------------------------*/
/* Attach file handle for energy file #energyfile_                      */
void STR::TimInt::AttachEnergyFile(std::string name)
{
  if (not energyfile_ or name != "")
  {
    // if energy file with new name is attached, delete the old file handle
    DetachEnergyFile();
    std::string energyname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()
      + name + ".energy";
    energyfile_ = new std::ofstream(energyname.c_str());
    *energyfile_ << "# timestep time total_energy"
                 << " kinetic_energy internal_energy external_energy"
                 << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*/
/* Return (rotatory) transformation matrix of local co-ordinate systems */
Teuchos::RCP<const LINALG::SparseMatrix> STR::TimInt::GetLocSysTrafo() const
{
  if (locsysman_ != Teuchos::null)
    return locsysman_->Trafo();

  return Teuchos::null;
}

//! Return stiffness,
//! i.e. force residual differentiated by displacements
Teuchos::RCP<LINALG::SparseMatrix> STR::TimInt::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_);
}

//! Return stiffness,
//! i.e. force residual differentiated by displacements
Teuchos::RCP<LINALG::BlockSparseMatrixBase> STR::TimInt::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(stiff_);
}

//! Return sparse mass matrix
Teuchos::RCP<LINALG::SparseMatrix> STR::TimInt::MassMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mass_);
}


/*----------------------------------------------------------------------*/
/* Return domain map of mass matrix                                     */
const Epetra_Map& STR::TimInt::GetDomainMap()
{
  return mass_->DomainMap();
}

/*----------------------------------------------------------------------*/
/* Creates the field test                                               */
Teuchos::RCP<DRT::ResultTest> STR::TimInt::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(*this));
}

/*----------------------------------------------------------------------*/
/* dof map of vector of unknowns                                        */
Teuchos::RCP<const Epetra_Map> STR::TimInt::DofRowMap()
{
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* dof map of vector of unknowns                                        */
/* new method for multiple dofsets                                      */
Teuchos::RCP<const Epetra_Map> STR::TimInt::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------
 Shorten the Dirichlet DOF set
 ----------------------------------------------------------------------*/
void STR::TimInt::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > othermaps;
  othermaps.push_back(maptoremove);
  othermaps.push_back(dbcmaps_->OtherMap());
  Teuchos::RCP<Epetra_Map> othermerged = LINALG::MultiMapExtractor::MergeMaps(othermaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), othermerged, false);
  return;
}

/*----------------------------------------------------------------------*/
/* reset everything (needed for biofilm simulations)                    */
void STR::TimInt::Reset()
{
  // displacements D_{n}
  dis_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // displacements D_{n}
  dism_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));

  // displacements D_{n+1} at t_{n+1}
  disn_ = LINALG::CreateVector(*dofrowmap_, true);
  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*dofrowmap_, true);
  // accelerations A_{n+1} at t_{n+1}
  accn_ = LINALG::CreateVector(*dofrowmap_, true);
  // create empty interface force vector
  fifc_ = LINALG::CreateVector(*dofrowmap_, true);

  // set initial fields
  SetInitialFields();

  return;
}

/*----------------------------------------------------------------------*/
/* set structure displacement vector due to biofilm growth              */
void STR::TimInt::SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp)
{
  strgrdisp_= struct_growth_disp;

  return;
}

/*----------------------------------------------------------------------*/
/* Resize MStep Object due to time adaptivity in FSI                    */
void STR::TimInt::ResizeMStepTimAda()
{
  // resize time and stepsize fields
  time_->Resize(-1, 0, (*time_)[0]);
  dt_->Resize(-1, 0, (*dt_)[0]);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  dis_->Resize(-1, 0, dofrowmap_, true);
  vel_->Resize(-1, 0, dofrowmap_, true);
  acc_->Resize(-1, 0, dofrowmap_, true);

  return;
}
