/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_implicit.cpp
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme

     o generalized-alpha time-integration scheme

     o implicit characteristic Galerkin (ICG) time-integration scheme (only for level-set transport)

     o explicit taylor galerkin (TG) time-integration schemes (only for level-set transport)

     and stationary solver.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_implicit.H"
#include "scatra_ele_action.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid/fluid_utils.H" // for splitter
#include "scatra_utils.H"
#include "scatra_utils_splitstrategy.H" // for blockmatrix-splitstrategy
#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"
#include "../drt_fluid/dyn_smag.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_SerialDenseVector.h>

#include "../drt_fluid/fluid_meshtying.H"

// for the condition writer output
/*
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
*/
/*
// for output of intermediate states in Newton loop
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
*/

//#define VISUALIZE_ELEDATA_GMSH
//only if VISUALIZE_ELEDATA_GMSH
//#include "../drt_io/io_gmsh.H"


/*==========================================================================*/
// Constructors and destructors and related methods
/*==========================================================================*/

/*----------------------------------------------------------------------*
 |  Constructor                                        (public) vg 05/07|
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::ScaTraTimIntImpl(
    Teuchos::RCP<DRT::Discretization>      actdis,
    Teuchos::RCP<LINALG::Solver>           solver,
    Teuchos::RCP<Teuchos::ParameterList>   params,
    Teuchos::RCP<Teuchos::ParameterList>   extraparams,
    Teuchos::RCP<IO::DiscretizationWriter> output) :
  // call constructor for "nontrivial" objects
  solver_ (solver),
  params_ (params),
  extraparams_(extraparams),
  myrank_ (actdis->Comm().MyPID()),
  // splitter, // not initialized
  errfile_  (extraparams->get<FILE*>("err file")),
  scatratype_  (DRT::INPUT::IntegralValue<INPAR::SCATRA::ScaTraType>(*params,"SCATRATYPE")),
  isale_    (extraparams->get<bool>("isale")),
  solvtype_ (DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params,"SOLVERTYPE")),
  // incremental_, // not initialized
  project_(false),
  initialvelset_(false),
  fssgd_    (DRT::INPUT::IntegralValue<INPAR::SCATRA::FSSUGRDIFF>(*params,"FSSUGRDIFF")),
  // turbmodel_, // not initialized
  writeflux_(DRT::INPUT::IntegralValue<INPAR::SCATRA::FluxType>(*params,"WRITEFLUX")),
  // writefluxids_, // not initialized
  outmean_  (DRT::INPUT::IntegralValue<int>(*params,"OUTMEAN")),
  outputgmsh_(DRT::INPUT::IntegralValue<int>(*params,"OUTPUT_GMSH")),
  time_   (0.0),
  maxtime_  (params->get<double>("MAXTIME")),
  step_   (0),
  stepmax_  (params->get<int>("NUMSTEP")),
  dta_      (params->get<double>("TIMESTEP")),
  dtele_(0.0),
  dtsolve_(0.0),
  timealgo_ (DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(*params,"TIMEINTEGR")),
  numscal_(0),
  // phi vectors, // all not initialized
  // vel_, // not initialized
  // many many not initialized
  cdvel_    (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*params,"VELOCITYFIELD")),
  discret_(actdis),
  output_ (output),
  convform_ (DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(*params,"CONVFORM")),
  // many many not initialized
  lastfluxoutputstep_(-1),
  msht_(DRT::INPUT::IntegralValue<int>(*params,"MESHTYING")),
  gstatnumite_(0),
  gstatincrement_(0.0),
  frt_      (0.0),
  numinflowsteps_(extraparams->sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")),
  turbinflow_(DRT::INPUT::IntegralValue<int>(extraparams->sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")),
  reinitswitch_(extraparams->get<bool>("REINITSWITCH",false)),
  w_(Teuchos::null),
  c_(Teuchos::null),
  upres_    (params->get<int>("UPRES")),
  uprestart_(params->get<int>("RESTARTEVRY")),
  neumanninflow_(DRT::INPUT::IntegralValue<int>(*params,"NEUMANNINFLOW")),
  convheatrans_(DRT::INPUT::IntegralValue<int>(*params,"CONV_HEAT_TRANS")),
  skipinitder_(DRT::INPUT::IntegralValue<int>(*params,"SKIPINITDER"))
{
  // what kind of equations do we actually want to solve?
  // (For the moment, we directly conclude from the problem type, Only ELCH applications
  //  allow the usage of a given user input)
  // additional exception: turbulent passive scalar transport: only for this case and loma
  // vectors and variables for turbulence models are provided
  PROBLEM_TYP prbtype = DRT::Problem::Instance()->ProblemType();
  if ((scatratype_ == INPAR::SCATRA::scatratype_undefined) or
     ((prbtype != prb_elch) and (scatratype_ != INPAR::SCATRA::scatratype_turbpassivesca)))
  {
    if (prbtype == prb_elch)              scatratype_ = INPAR::SCATRA::scatratype_elch_enc;
    else if (prbtype == prb_combust)      scatratype_ = INPAR::SCATRA::scatratype_levelset;
    else if (prbtype == prb_loma)         scatratype_ = INPAR::SCATRA::scatratype_loma;
    else if (prbtype == prb_scatra)       scatratype_ = INPAR::SCATRA::scatratype_condif;
    else if (prbtype == prb_gas_fsi)      scatratype_ = INPAR::SCATRA::scatratype_condif;
    else if (prbtype == prb_biofilm_fsi)  scatratype_ = INPAR::SCATRA::scatratype_condif;
    else if (prbtype == prb_thermo_fsi)   scatratype_ = INPAR::SCATRA::scatratype_loma;
    else if (prbtype == prb_poroscatra)   scatratype_ = INPAR::SCATRA::scatratype_poro;
    else if (prbtype == prb_ssi)          scatratype_ = INPAR::SCATRA::scatratype_condif;
    else
      dserror("Problemtype %s not supported", DRT::Problem::Instance()->ProblemName().c_str());
  }

  // -------------------------------------------------------------------
  // determine whether linear incremental or nonlinear solver
  // -------------------------------------------------------------------
  switch(solvtype_)
  {
  case INPAR::SCATRA::solvertype_nonlinear:
  case INPAR::SCATRA::solvertype_linear_incremental:
  {
    incremental_ = true;
  }
  break;
  case INPAR::SCATRA::solvertype_linear_full:
  {
    incremental_ = false;
  }
  break;
  default:
    dserror("Received illegal scatra solvertype enum.");
    break;
  }

  // -------------------------------------------------------------------
  // check compatibility of boundary conditions
  // -------------------------------------------------------------------
  if (neumanninflow_ and convheatrans_)
    dserror("Neumann inflow and convective heat transfer boundary conditions must not appear simultaneously for the same problem!");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  // do not apply periodic boundary conditions for the second scatra reinit object
  // then again a new redistribution of the redistributed scatra discretization would be performed
  if(reinitswitch_ == false)
  {
    pbc_ = Teuchos::rcp(new PeriodicBoundaryConditions (discret_, false));
    pbc_->UpdateDofsForPeriodicBoundaryConditions();
    pbcmapmastertoslave_ = pbc_->ReturnAllCoupledRowNodes();
  }


  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs()))
    discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  int mynumscal= 0;
  if (discret_->NumMyRowNodes()>0)
    mynumscal = discret_->NumDof(0,discret_->lRowNode(0));
  // to support completely empty procs, communication is required
  discret_->Comm().MaxAll(&mynumscal,&numscal_,1);

  if (IsElch(scatratype_))
  {
    if (numscal_ > 1) // we have at least two ion species + el. potential
    {
      // number of concentrations transported is numdof-1
      numscal_ -= 1;

      Teuchos::ParameterList& diffcondparams = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

      // currrent is a solution variable
      if(DRT::INPUT::IntegralValue<int>(diffcondparams,"CURRENT_SOLUTION_VAR"))
      {
        // shape of local row element(0) -> number of space dimensions
        //int dim = DRT::Problem::Instance()->NDim();
        int dim = DRT::UTILS::getDimension(discret_->lRowElement(0)->Shape());
        // number of concentrations transported is numdof-1-nsd
        numscal_ -= dim;
      }

      if(DRT::INPUT::IntegralValue<int>(diffcondparams,"DIFFCOND_FORMULATION"))
        ValidParameterDiffCond();

    }
    // set up the concentration-el.potential splitter
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_,*splitter_);
  }
  else if (scatratype_ == INPAR::SCATRA::scatratype_loma and numscal_ > 1)
  {
    // set up a species-temperature splitter (if more than one scalar)
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_-1,*splitter_);
  }

  if (scatratype_ == INPAR::SCATRA::scatratype_turbpassivesca and numscal_ > 1)
   dserror("Turbulent passive scalar transport not supported for more than one scalar!");

  if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND")
      and msht_ == INPAR::FLUID::no_meshtying)
  {
    // we need a block sparse matrix here
    if (not IsElch(scatratype_))
      dserror("Block-Preconditioning is only for ELCH problems");

    if ((scatratype_!=INPAR::SCATRA::scatratype_elch_enc))
      dserror("Special ELCH assemble strategy for block-matrix will not assemble A_11 block!");

    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal_ due to electroneutrality, A_11: EMPTY matrix !!!!!
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    Teuchos::RCP<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(*splitter_,*splitter_,27,false,true));
    blocksysmat->SetNumScal(numscal_);

    sysmat_ = blocksysmat;
  }
  else if(msht_!= INPAR::FLUID::no_meshtying)
  {
    if (((msht_== INPAR::FLUID::condensed_bmat) or
        (msht_== INPAR::FLUID::condensed_bmat_merged) or
         msht_== INPAR::FLUID::coupling_iontransport_laplace) and
        (scatratype_ == INPAR::SCATRA::scatratype_elch_enc))
      dserror("In the context of mesh-tying, the ion-transport system inluding the electroneutrality condition \n"
          "cannot be solved by a block matrix");

    if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND") and
        (msht_== INPAR::FLUID::condensed_bmat or
                msht_== INPAR::FLUID::condensed_bmat_merged))
      dserror("Switch of the block matrix!!");

    // define parameter list for meshtying
    Teuchos::ParameterList mshtparams;
    mshtparams.set("theta", params_->get<double>("THETA"));
    mshtparams.set<int>("mshtoption", msht_);

    meshtying_ = Teuchos::rcp(new FLD::Meshtying(discret_, *solver_, mshtparams));
    sysmat_ = meshtying_->Setup();
  }
  else
  {
    // initialize standard (stabilized) system matrix (and save its graph!)
    // in standard case, but do not save the graph if fine-scale subgrid
    // diffusivity is used in non-incremental case
    if (fssgd_ != INPAR::SCATRA::fssugrdiff_no and not incremental_)
    {
      // this is a very special case
      // only fssugrdiff_artificial is allowed in combination with non-incremental
      sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));
    }
    else sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
  }

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1 and n
  phinp_ = LINALG::CreateVector(*dofrowmap,true);
  phin_  = LINALG::CreateVector(*dofrowmap,true);

  if(reinitswitch_)
  {
    phistart_  = LINALG::CreateVector(*dofrowmap,true);
  }

  // temporal solution derivative at time n+1
  phidtnp_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  hist_ = LINALG::CreateVector(*dofrowmap,true);

  // velocities (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  convel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
  vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // acceleration and pressure required for computation of subgrid-scale
  // velocity (always four components per node)
  accpre_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,4,true));

  if (isale_)
  {
    // displacement field for moving mesh applications using ALE
    // (get noderowmap of discretization for creating this multivector)
    dispnp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
  }

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  // -------------------------------------------------------------------
  // ensure that the Transport string was removed from conditions
  // -------------------------------------------------------------------
  {
    DRT::Condition* cond = discret_->GetCondition("TransportDirichlet");
    if (cond) dserror("Found a Transport Dirichlet condition. Remove Transport string!");
    cond = discret_->GetCondition("TransportNeumann");
    if (cond) dserror("Found a Transport Neumann condition. Remove Transport string!");
  }

  // -------------------------------------------------------------------
  // create vectors associated to solution process
  // -------------------------------------------------------------------
  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // the residual vector --- more or less the rhs
  residual_ = LINALG::CreateVector(*dofrowmap,true);

  // residual vector containing the normal boundary fluxes
  trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

  // incremental solution vector
  increment_ = LINALG::CreateVector(*dofrowmap,true);

  // subgrid-diffusivity(-scaling) vector
  // (used either for AVM3 approach or temperature equation
  //  with all-scale subgrid-diffusivity model)
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    subgrdiff_ = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // set parameters associated to potential statistical flux evaluations
  // -------------------------------------------------------------------
  // get fluid turbulence sublist
  Teuchos::ParameterList * turbparams =&(extraparams_->sublist("TURBULENCE MODEL"));

  // parameters for statistical evaluation of normal fluxes
  samstart_  = turbparams->get<int>("SAMPLING_START");
  samstop_   = turbparams->get<int>("SAMPLING_STOP" );
  dumperiod_ = turbparams->get<int>("DUMPING_PERIOD");
  if (dumperiod_ < 0) dserror("dumperiod_ is negative!");

  // initialize vector for statistics (assume a maximum of 10 conditions)
  sumnormfluxintegral_ = Teuchos::rcp(new Epetra_SerialDenseVector(10));

  // get desired scalar id's for flux output
  {
    int    word1;
    std::istringstream mystream(Teuchos::getNumericStringParameter(*params_,"WRITEFLUX_IDS"));
    while (mystream >> word1)
      writefluxids_.push_back(word1);

    if (writefluxids_[0]==(-1)) //default is to perform flux output for ALL scalars
    {
      writefluxids_.resize(numscal_);
      for(int k=0;k<numscal_;++k)
        writefluxids_[k]=k+1;
    }

    if ((writeflux_!=INPAR::SCATRA::flux_no) and (myrank_ == 0))
    {
      IO::cout << "Flux output is performed for scalars: ";
      for (unsigned int i=0; i < writefluxids_.size();i++)
      {
        const int id = writefluxids_[i];
        IO::cout << writefluxids_[i] << " ";
        if ((id<1) or (id > numscal_)) // check validity of these numbers as well !
          dserror("Received illegal scalar id for flux output: %d",id);
      }
      IO::cout << IO::endl;
    }
  }

  // -------------------------------------------------------------------
  // necessary only for AVM3 approach:
  // initialize subgrid-diffusivity matrix + respective output
  // -------------------------------------------------------------------
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
  {
    sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

    // fine-scale velocities (always three velocity components per node)
    // transferred from the fluid field
    // only Smagorinsky small
    if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
      fsvel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));

    // Output
    if (myrank_ == 0)
    {
      cout << "Fine-scale subgrid-diffusivity approach based on AVM3: ";
      cout << fssgd_;
      cout << &endl << &endl;
    }

    if (scatratype_ == INPAR::SCATRA::scatratype_loma or
        scatratype_ == INPAR::SCATRA::scatratype_turbpassivesca)
    {
      if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small
          and turbparams->get<string>("FSSUGRVISC") != "Smagorinsky_small")
        dserror ("Same subgrid-viscosity approach expected!");
      if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all
          and turbparams->get<string>("FSSUGRVISC") != "Smagorinsky_all")
        dserror ("Same subgrid-viscosity approach expected!");
    }
  }

  // -------------------------------------------------------------------
  // get turbulence model and parameters for low-Mach-number case
  // -------------------------------------------------------------------
  turbmodel_ = INPAR::FLUID::no_model;
  if (scatratype_ == INPAR::SCATRA::scatratype_loma or
      scatratype_ == INPAR::SCATRA::scatratype_turbpassivesca)
  {
    // set turbulence model
    if (turbparams->get<string>("PHYSICAL_MODEL") == "Smagorinsky")
    {
      turbmodel_ = INPAR::FLUID::smagorinsky;

      // Output
      if (turbmodel_ and myrank_ == 0)
      {
        cout << "All-scale subgrid-diffusivity model: ";
        cout << turbparams->get<string>("PHYSICAL_MODEL");
        cout << &endl << &endl;
      }
    }
    else if (turbparams->get<string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
    {
      turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
      // access to the dynamic Smagorinsky class will provided by the
      // scatra fluid couling algorithm
    }
    else if (turbparams->get<string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
    {
      turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;

      // initalize matrix used to build the scale separation operator
      sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

      // fine-scale velocities (always three velocity components per node)
      // transferred from the fluid field
      fsvel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));

      Teuchos::ParameterList * mfsparams =&(extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES"));
      if (mfsparams->get<string>("SCALE_SEPARATION")!= "algebraic_multigrid_operator")
       dserror("Only scale separation by plain algebraic multigrid available in scatra!");

      // Output
      if (turbmodel_ and myrank_ == 0)
      {
        cout << "Multifractal subgrid-scale model: ";
        cout << turbparams->get<string>("PHYSICAL_MODEL");
        cout << &endl << &endl;
      }
    }

    // warning No. 1: if classical (all-scale) turbulence model other than
    // Smagorinsky or multifractal subrgid-scale modeling
    // is intended to be used
    if (turbparams->get<string>("PHYSICAL_MODEL") != "Smagorinsky" and
        turbparams->get<string>("PHYSICAL_MODEL") != "Dynamic_Smagorinsky" and
        turbparams->get<string>("PHYSICAL_MODEL") != "Multifractal_Subgrid_Scales" and
        turbparams->get<string>("PHYSICAL_MODEL") != "no_model")
      dserror("No classical (all-scale) turbulence model other than constant-coefficient Smagorinsky model and multifractal subrgid-scale modeling currently possible!");

    // warning No. 2: if classical (all-scale) turbulence model and fine-scale
    // subgrid-viscosity approach are intended to be used simultaneously
    if (turbmodel_==INPAR::FLUID::smagorinsky and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
      dserror("No combination of classical turbulence model and fine-scale subgrid-diffusivity approach currently possible!");
  }

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(*params_,"INITIALFIELD"),
      params_->get<int>("INITFUNCNO"));

  // initializes variables for natural convection (ELCH) if necessary
  SetupElchNatConv();

  // screen output (has to come after SetInitialField)
  if (IsElch(scatratype_))
  {
    // a safety check for the solver type
    if ((numscal_ > 1) && (solvtype_!=INPAR::SCATRA::solvertype_nonlinear))
      dserror("Solver type has to be set to >>nonlinear<< for ion transport.");

    frt_ = INPAR::SCATRA::faraday_const/(8.314472 * extraparams_->get<double>("TEMPERATURE"));

    if (myrank_==0)
    {
      cout<<"\nSetup of splitter: numscal = "<<numscal_<<endl;
      cout<<"Temperature value T (Kelvin)     = "<<extraparams_->get<double>("TEMPERATURE")<<endl;
      cout<<"Constant F/RT                    = "<<frt_<<endl;
    }

    // setup of magnetic field (always three(!) components per node)
    const int magnetfuncno = (extraparams_->sublist("ELCH CONTROL")).get<int>("MAGNETICFIELD_FUNCNO");
    if (magnetfuncno > 0)
    {
      // allocate the multivector
      magneticfield_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
      // fill it with values
      SetMagneticField(magnetfuncno);
    }
  }
  if(extraparams_->isSublist("ELCH CONTROL")) // did we get the elch sublist?
  {
    // conductivity must be stored for the galvanostatic condition in a global variable
    sigma_ = ComputeConductivity(); // every processor has to do this call
    if (myrank_==0)
    {
      for (int k=0;k < numscal_;k++)
      {
        cout<<"Electrolyte conductivity (species "<<k+1<<")    = "<<sigma_[k]<<endl;
      }
      if (scatratype_==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
      {
        double diff = sigma_[0];
        for (int k=1;k < numscal_;k++)
        {
          diff += sigma_[k];
        }
        cout<<"Electrolyte conductivity (species elim) = "<<sigma_[numscal_]-diff<<endl;
      }
      cout<<"Electrolyte conductivity (all species)  = "<<sigma_[numscal_]<<endl<<endl;
    }
  }

  // sysmat might be singular (some modes are defined only up to a constant)
  // in this case, we need basis vectors for the nullspace/kernel
  vector<DRT::Condition*> KSPCond;
  discret_->GetCondition("KrylovSpaceProjection",KSPCond);
  int numcond = KSPCond.size();
  int nummodes = 0;
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPCond[icond]->Get<std::string>("discretization");
    if (*name == "scatra") nummodes++;
  }
  if (nummodes > 0)
  {
    project_ = true;
    PrepareKrylovSpaceProjection();
    if (myrank_ == 0)
      cout<<"\nSetup of KrylovSpaceProjection:\n"<<
      "    => number of kernel basis vectors: "<<nummodes<<endl<<endl;
  }

  return;

} // ScaTraTimIntImpl::ScaTraTimIntImpl

/*----------------------------------------------------------------------*
 | Destructor dtor                                   (public) gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl()
{
  return;
}

/*==========================================================================*/
// general framework
/*==========================================================================*/

/*--- set, prepare, and predict --------------------------------------------*/

/*----------------------------------------------------------------------*
 | setup the variables to do a new time step          (public)  vg 08/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareTimeStep()
{
  // time measurement: prepare time step
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + prepare time step");

  // -------------------------------------------------------------------
  //                       initialization
  // -------------------------------------------------------------------
  if (step_ == 0 and not skipinitder_)
  {
    // if initial velocity field has not been set here, the initial time derivative of phi will be
    // calculated wrongly for some time integration schemes
    if (initialvelset_) PrepareFirstTimeStep();
    else if (reinitswitch_){}
    else dserror("Initial velocity field has not been set");
  }

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  IncrementTimeAndStep();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  ApplyDirichletBC(time_,phinp_,Teuchos::null);
  ApplyNeumannBC(time_,phinp_,neumann_loads_);

  // -------------------------------------------------------------------
  //     update velocity field if given by function AND time curve
  // -------------------------------------------------------------------
  if (cdvel_ == INPAR::SCATRA::velocity_function_and_curve)
    SetVelocityField();

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if ((step_==1 or (turbinflow_ and step_==numinflowsteps_+1)) and
      (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
     AVM3Preparation();

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  ComputeIntermediateValues();

  return;

} // ScaTraTimIntImpl::PrepareTimeStep

/*----------------------------------------------------------------------*
 | preparations for solve                                (public) mr.x  |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareLinearSolve()
{
  // special preparations for multifractal subgrid-scale model
  if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales)
    RecomputeMeanCsgsB();
  
  // call elements to calculate system matrix and rhs and assemble
  AssembleMatAndRHS();

  // potential residual scaling and potential addition of Neumann terms
  ScalingAndNeumann();

  // apply Dirichlet boundary conditions
  ApplyDirichletToSystem();
}

/*----------------------------------------------------------------------*
 | update the velocity field                                  gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField()
{
  if (cdvel_ == INPAR::SCATRA::velocity_zero)
  {
    convel_->PutScalar(0.); // just to be sure!
    vel_->PutScalar(0.);
  }
  else if ((cdvel_ == INPAR::SCATRA::velocity_function)
      or (cdvel_ == INPAR::SCATRA::velocity_function_and_curve))
  {
    int err(0);
    const int numdim = 3; // the velocity field is always 3D
    const int velfuncno = params_->get<int>("VELFUNCNO");
    const int velcurveno = params_->get<int>("VELCURVENO");
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value = DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(index,lnode->X(),time_,NULL);
        if (cdvel_ == INPAR::SCATRA::velocity_function_and_curve)
        {
          value *= DRT::Problem::Instance()->Curve(velcurveno-1).f(time_);
        }
        err = convel_->ReplaceMyValue (lnodeid, index, value);
        if (err!=0) dserror("error while inserting a value into convel_");
        err = vel_->ReplaceMyValue (lnodeid, index, value);
        if (err!=0) dserror("error while inserting a value into vel_");
      }
    }
  }
  else
    dserror("Wrong SetVelocity() action for velocity field type %d!",cdvel_);

  // initial velocity field has now been set
  if (step_ == 0) initialvelset_ = true;

  return;

} // ScaTraImplicitTimeInt::SetVelocityField

/*----------------------------------------------------------------------*
 | set convective velocity field (+ pressure and acceleration field as  |
 | well as fine-scale velocity field, if required)            gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(
Teuchos::RCP<const Epetra_Vector> convvel,
Teuchos::RCP<const Epetra_Vector> acc,
Teuchos::RCP<const Epetra_Vector> vel,
Teuchos::RCP<const Epetra_Vector> fsvel,
Teuchos::RCP<const DRT::DofSet>   dofset,
Teuchos::RCP<DRT::Discretization> dis)
{
  //---------------------------------------------------------------------------
  // preliminaries
  //---------------------------------------------------------------------------
  if (cdvel_ != INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Wrong SetVelocityField() called for velocity field type %d!",cdvel_);

  TEUCHOS_FUNC_TIME_MONITOR("SCATRA: set convective velocity field");

//#ifdef DEBUG   // is this costly, when we do this test always?
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not dis->NodeRowMap()->SameAs(*(discret_->NodeRowMap())))
    dserror("Fluid/Structure and Scatra noderowmaps are NOT identical. Emergency!");
//#endif

  // define error variable
  int err(0);

  // boolean indicating whether acceleration vector exists
  // -> if yes, subgrid-scale velocity may need to be computed on element level
  bool sgvelswitch = (acc != Teuchos::null);

  // boolean indicating whether fine-scale velocity vector exists
  // -> if yes, multifractal subgrid-scale modeling is applied
  bool fsvelswitch = (fsvel != Teuchos::null);

  // some thing went wrong if we want to use multifractal subgrid-scale modeling
  // and have not got the fine-scale velocity
  if (step_>=1 and (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales
       or fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
       and not fsvelswitch)
    dserror("Fine-scale velocity expected for multifractal subgrid-scale modeling!");
  // as fsvelswitch is also true for smagorinsky_all, we have to reset fsvelswitch
  // as the corresponding vector, which is not necessary, is not provided in scatra
  if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all and fsvelswitch)
    fsvelswitch = false;

  //---------------------------------------------------------------------------
  // transfer of dofs
  // (We rely on the fact that the scatra discretization is a clone of the
  // fluid or structure mesh, respectively, meaning that a scatra node has the
  // same local (and global) ID as its corresponding fluid/structure node.)
  //---------------------------------------------------------------------------
  // loop over all local nodes of scatra discretization
  for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // get local fluid/structure node with the same lnodeid
    DRT::Node* lnode = dis->lRowNode(lnodeid);

    // care for the slave nodes of rotationally symm. periodic boundary conditions
    double rotangle(0.0);
    bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(lnode,rotangle);

    // get degrees of freedom associated with this fluid/structure node
    // two particular cases have to be considered:
    // - in non-XFEM case, the first dofset is always considered, allowing for
    //   using multiple dof sets, e.g., for structure-based scalar transport
    // - for XFEM, a different nodeset is required
    vector<int> nodedofs;
    if (dofset == Teuchos::null) nodedofs = dis->Dof(0,lnode);
    else                         nodedofs = (*dofset).Dof(lnode);

    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    //-------------------------------------------------------------------------
    // transfer of velocity dofs
    //-------------------------------------------------------------------------
    for (int index=0;index < numdim; ++index)
    {
      // get global and local ID
      const int gid = nodedofs[index];
      // const int lid = dofrowmap->LID(gid);
      const int lid = convvel->Map().LID(gid);
      if (lid < 0) dserror("Local ID not found in map for given global ID!");

      //-----------------------------------------------------------------------
      // get convective velocity
      //-----------------------------------------------------------------------
      double convelocity = (*convvel)[lid];

      // component of rotated vector field
      if (havetorotate)  convelocity = FLD::GetComponentOfRotatedVectorField(index,convvel,lid,rotangle);

      // insert velocity value into node-based vector
      err = convel_->ReplaceMyValue(lnodeid,index,convelocity);
      if (err != 0) dserror("Error while inserting value into vector convel_!");

      //-----------------------------------------------------------------------
      // get velocity
      //-----------------------------------------------------------------------
      if (vel != Teuchos::null)
      {
        // get value of corresponding velocity component
        double velocity = (*vel)[lid];

        // component of rotated vector field
        if (havetorotate) velocity = FLD::GetComponentOfRotatedVectorField(index,vel,lid,rotangle);

        // insert velocity value into node-based vector
        err = vel_->ReplaceMyValue(lnodeid,index,velocity);
        if (err != 0) dserror("Error while inserting value into vector vel_!");
      }
      else
      {
        // if velocity vector is not provided by the respective algorithm, we
        // assume that it equals the given convective velocity:
        // insert velocity value into node-based vector
        err = vel_->ReplaceMyValue(lnodeid,index,convelocity);
        if (err != 0) dserror("Error while inserting value into vector vel_!");
      }

      //-----------------------------------------------------------------------
      // get acceleration, if required
      //-----------------------------------------------------------------------
      if (sgvelswitch)
      {
        // get value of corresponding acceleration component
        double acceleration = (*acc)[lid];

        // component of rotated vector field
        if (havetorotate) acceleration = FLD::GetComponentOfRotatedVectorField(index,acc,lid,rotangle);

        // insert acceleration value into node-based vector
        err = accpre_->ReplaceMyValue(lnodeid,index,acceleration);
        if (err != 0) dserror("Error while inserting value into vector accpre_!");
      }

      //-----------------------------------------------------------------------
      // get fine-scale velocity, if required
      //-----------------------------------------------------------------------
      if (fsvelswitch)
      {
        // get value of corresponding fine-scale velocity component
        double fsvelocity = (*fsvel)[lid];

        // component of rotated vector field
        if (havetorotate) fsvelocity = FLD::GetComponentOfRotatedVectorField(index,fsvel,lid,rotangle);

        // insert fine-scale velocity value into node-based vector
        err = fsvel_->ReplaceMyValue(lnodeid,index,fsvelocity);
        if (err != 0) dserror("Error while inserting value into vector fsvel_!");
      }
    }

    //-------------------------------------------------------------------------
    // transfer of pressure dofs, if required
    //-------------------------------------------------------------------------
    if (sgvelswitch)
    {
      // get global and local ID
      const int gid = nodedofs[numdim];
      // const int lid = dofrowmap->LID(gid);
      const int lid = convvel->Map().LID(gid);
      if (lid < 0) dserror("Local ID not found in map for given global ID!");

      // get value of corresponding pressure component
      double pressure = (*convvel)[lid];

      // insert pressure value into node-based vector
      err = accpre_->ReplaceMyValue(lnodeid,numdim,pressure);
      if (err != 0) dserror("Error while inserting value into vector accpre_!");
    }

    //-------------------------------------------------------------------------
    // to be sure for 1- and 2-D problems:
    // set all unused velocity components to zero
    //-------------------------------------------------------------------------
    for (int index=numdim; index < 3; ++index)
    {
      err = convel_->ReplaceMyValue(lnodeid,index,0.0);
      if (err != 0) dserror("Error while inserting value into vector convel_!");

      err = vel_->ReplaceMyValue(lnodeid,index,0.0);
      if (err != 0) dserror("Error while inserting value into vector vel_!");
    }
  }

  // confirm that initial velocity field has now been set
  if (step_ == 0) initialvelset_ = true;

  return;

} // ScaTraTimIntImpl::SetVelocityField

/*----------------------------------------------------------------------------*
 | Redistribute the scatra discretization and vectors         rasthofer 07/11 |
 | according to nodegraph according to nodegraph              DA wichmann     |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  if (reinitswitch_ == false)
  {
    const Epetra_BlockMap oldphinpmap = phinp_->Map();

    // the rowmap will become the new distribution of nodes
    const Epetra_BlockMap rntmp = nodegraph->RowMap();
    Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,discret_->Comm());

    // the column map will become the new ghosted distribution of nodes
    const Epetra_BlockMap Mcntmp = nodegraph->ColMap();
    Epetra_Map newnodecolmap(-1,Mcntmp.NumMyElements(),Mcntmp.MyGlobalElements(),0,discret_->Comm());

    // do the redistribution
    discret_->Redistribute(newnoderowmap,newnodecolmap, false, false, false);

    // update the PBCs and PBCDofSet
    pbc_->PutAllSlavesToMastersProc();
    pbcmapmastertoslave_ = pbc_->ReturnAllCoupledRowNodes();

    // ensure that degrees of freedom in the discretization have been set
    if ((not discret_->Filled()) or (not discret_->HaveDofs()))
      discret_->FillComplete();
  }

  //--------------------------------------------------------------------
  // Now update all Epetra_Vectors and Epetra_Matrix to the new dofmap
  //--------------------------------------------------------------------

  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  if (IsElch(scatratype_))
  {
    // set up the concentration-el.potential splitter
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_,*splitter_);
  }
  else if (scatratype_ == INPAR::SCATRA::scatratype_loma and numscal_ > 1)
  {
    // set up a species-temperature splitter (if more than one scalar)
    splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
    FLD::UTILS::SetupFluidSplit(*discret_,numscal_-1,*splitter_);
  }

  if (scatratype_ == INPAR::SCATRA::scatratype_turbpassivesca and numscal_ > 1)
   dserror("Turbulent passive scalar transport not supported for more than one scalar!");

  if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND")
      and msht_ == INPAR::FLUID::no_meshtying)
  {
    // we need a block sparse matrix here
    if (not IsElch(scatratype_))
      dserror("Block-Preconditioning is only for ELCH problems");
    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal_ due to electroneutrality, A_11: empty matrix
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    Teuchos::RCP<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(*splitter_,*splitter_,27,false,true));
    blocksysmat->SetNumScal(numscal_);

    sysmat_ = blocksysmat;
  }
  else if(msht_!= INPAR::FLUID::no_meshtying)
  {
    if (msht_!= INPAR::FLUID::condensed_smat)
      dserror("In the moment the only option is condensation in a sparse matrix");

    // define parameter list for meshtying
    Teuchos::ParameterList mshtparams;
    mshtparams.set("theta", params_->get<double>("THETA"));
    mshtparams.set<int>("mshtoption", msht_);

    meshtying_ = Teuchos::rcp(new FLD::Meshtying(discret_, *solver_, mshtparams));
    sysmat_ = meshtying_->Setup();
  }
  else
  {
    // initialize standard (stabilized) system matrix (and save its graph!)
    // in standard case, but do not save the graph if fine-scale subgrid
    // diffusivity is used in non-incremental case
    if (fssgd_ != INPAR::SCATRA::fssugrdiff_no and not incremental_)
         //sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));
         // cf constructor
         sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
    else sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));
  }

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------

  // solutions at time n+1 and n
  Teuchos::RCP<Epetra_Vector> old;

  if (phinp_ != Teuchos::null)
  {
    old = phinp_;
    phinp_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *phinp_);
  }

  if (phin_ != Teuchos::null)
  {
    old = phin_;
    phin_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *phin_);
  }

  // phi at time 0 as reference for reinitialization procedure
  if (phistart_ != Teuchos::null)
  {
    old = phistart_;
    phistart_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *phistart_);
  }

  // temporal solution derivative at time n+1
  if (phidtnp_ != Teuchos::null)
  {
    old = phidtnp_;
    phidtnp_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *phidtnp_);
  }

  // temporal solution derivative at time n
  if (phidtn_ != Teuchos::null)
  {
    old = phidtn_;
    phidtn_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *phidtn_);
  }

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  if (hist_ != Teuchos::null)
  {
    old = hist_;
    hist_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *hist_);
  }

  // velocities (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  Teuchos::RCP<Epetra_MultiVector> oldMulti;
  if (convel_ != Teuchos::null)
  {
    oldMulti = convel_;
    convel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
    LINALG::Export(*oldMulti, *convel_);
  }
  if (vel_ != Teuchos::null)
  {
    oldMulti = vel_;
    vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
    LINALG::Export(*oldMulti, *vel_);
  }
  if (fsvel_ != Teuchos::null)
  {
    oldMulti = fsvel_;
    fsvel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
    LINALG::Export(*oldMulti, *fsvel_);
  }

  // acceleration and pressure required for computation of subgrid-scale
  // velocity (always four components per node)
  if (accpre_ != Teuchos::null)
  {
    oldMulti = accpre_;
    accpre_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,4,true));
    LINALG::Export(*oldMulti, *accpre_);
  }

  if (dispnp_ != Teuchos::null)
  {
    oldMulti = dispnp_;
    dispnp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
    LINALG::Export(*oldMulti, *dispnp_);
  }

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  if (zeros_ != Teuchos::null)
  {
    old = zeros_;
    zeros_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *zeros_);
  }

  // -------------------------------------------------------------------
  // create vectors associated to solution process
  // -------------------------------------------------------------------
  // the vector containing body and surface forces
  if (neumann_loads_ != Teuchos::null)
  {
    old = neumann_loads_;
    neumann_loads_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *neumann_loads_);
  }

  // the residual vector --- more or less the rhs
  if (residual_ != Teuchos::null)
  {
    old = residual_;
    residual_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *residual_);
  }

  // residual vector containing the normal boundary fluxes
  if (trueresidual_ != Teuchos::null)
  {
    old = trueresidual_;
    trueresidual_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *trueresidual_);
  }

  // incremental solution vector
  if (increment_ != Teuchos::null)
  {
    old = increment_;
    increment_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *increment_);
  }

  // subgrid-diffusivity(-scaling) vector
  // (used either for AVM3 approach or temperature equation
  //  with all-scale subgrid-diffusivity model)
  if (subgrdiff_ != Teuchos::null)
  {
    old = subgrdiff_;
    subgrdiff_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *subgrdiff_);
  }

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
  {
    dserror("No redistribution for AVM3 subgrid stuff.");
  }

  if (IsElch(scatratype_))
    dserror("No redistribution for the elch.");

  if(discret_->Comm().MyPID()==0)
    cout << "done" << endl;

  return;
} // SCATRA::ScaTraTimIntImpl::Redistribute

/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::TimeLoop()
{
  // write out initial state
  // Output();

  // provide information about initial field (do not do for restarts!)
  if (Step()==0)
  {
    OutputElectrodeInfo();
    OutputMeanScalars();

    // compute error for problems with analytical solution (initial field!)
    EvaluateErrorComparedToAnalyticalSol();
  }

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  while ((step_<stepmax_) and ((time_+ EPS12) < maxtime_))
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                  solve nonlinear / linear equation
    // -------------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();

    // -------------------------------------------------------------------
    // evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  } // while

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();

  return;
} // ScaTraTimIntImpl::TimeLoop

/*----------------------------------------------------------------------*
 | contains the call of linear/nonlinear solver               gjb 02/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Solve()
{
  if (solvtype_==INPAR::SCATRA::solvertype_nonlinear)
    NonlinearSolve();
  else
    LinearSolve();
  //that's all
  return;
}

/*----------------------------------------------------------------------*
 | apply moving mesh data                                     gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> dispnp,
    Teuchos::RCP<DRT::Discretization> dis
)
{
  //---------------------------------------------------------------------------
  // only required in ALE case
  //---------------------------------------------------------------------------
  if (isale_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA: apply mesh movement");

    // check existence of displacement vector
    if (dispnp == Teuchos::null) dserror("Got null pointer for displacements!");

    // define error variable
    int err(0);

    // get dofrowmap of discretization
    const Epetra_Map* dofrowmap = dis->DofRowMap();

    //-------------------------------------------------------------------------
    // transfer of dofs
    // (We rely on the fact that the scatra discretization is a clone of the
    // fluid or structure mesh, respectively, meaning that a scatra node has the
    // same local (and global) ID as its corresponding fluid/structure node.)
    //-------------------------------------------------------------------------
    // loop over all local nodes of scatra discretization
    for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      // get local fluid/structure node with the same lnodeid
      DRT::Node* lnode = dis->lRowNode(lnodeid);

      // get degrees of freedom associated with this fluid/structure node
      // (first dofset always considered, allowing for using multiple
      //  dof sets, e.g., for structure-based scalar transport)
      vector<int> nodedofs = dis->Dof(0,lnode);

      // determine number of space dimensions
      const int numdim = DRT::Problem::Instance()->NDim();

      for (int index=0;index < numdim; ++index)
      {
        // get global and local ID
        const int gid = nodedofs[index];
        const int lid = dofrowmap->LID(gid);

        //---------------------------------------------------------------------
        // get displacement
        //---------------------------------------------------------------------
        double disp = (*dispnp)[lid];

        // insert displacement value into node-based vector
        err = dispnp_->ReplaceMyValue(lnodeid,index,disp);
        if (err != 0) dserror("Error while inserting value into vector dispnp_!");
      }

      //-----------------------------------------------------------------------
      // to be sure for 1- and 2-D problems:
      // set all unused displacement components to zero
      //-----------------------------------------------------------------------
      for (int index=numdim; index < 3; ++index)
      {
        err = dispnp_->ReplaceMyValue(lnodeid,index,0.0);
        if (err != 0) dserror("Error while inserting value into vector dispnp_!");
      }
    } // for lnodeid
  } // if (isale_)

  return;

} // ScaTraTimIntImpl::ApplyMeshMovement

/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector                        gjb   04/08|
 *----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFlux
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector field in comp. domain    gjb 06/09|
 *----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFluxInDomain
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 |  calculate mass / heat normal flux at specified boundaries  gjb 06/09|
 *----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFluxAtBoundary
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 |  print information about current time step to screen        mr. x    |
 *----------------------------------------------------------------------*/
inline void SCATRA::ScaTraTimIntImpl::PrintTimeStepInfo()
{
  if (myrank_==0)
    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d \n",
           time_,maxtime_,dta_,MethodTitle().c_str(),step_,stepmax_);
} // SCATRA::ScaTraTimIntImpl::PrintTimeStepInfo

/*----------------------------------------------------------------------*
 | return system matrix downcasted as sparse matrix           gjb 02/11 |
 | implemented here to be able to use forward declaration in .H         |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> SCATRA::ScaTraTimIntImpl::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------*
 | return system matrix downcasted as block sparse matrix     gjb 06/10 |
 | implemented here to be able to use forward declaration in .H         |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> SCATRA::ScaTraTimIntImpl::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                          gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Output()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // solution output and potentially restart data and/or flux data
  if (DoOutput())
  {
    // step number and time (only after that data output is possible)
    output_->NewStep(step_,time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_==upres_) output_->WriteElementData();

    // write state vectors
    OutputState();

    // write output to Gmsh postprocessing files
    if (outputgmsh_) OutputToGmsh(step_, time_);

    // add restart data
    if (step_%uprestart_==0) OutputRestart();

    // write flux vector field (only writing, calculation was done during Update() call)
    if (writeflux_!=INPAR::SCATRA::flux_no)
    {
      // for flux output of initial field (before first solve) do:
      if (step_==0)
        flux_=CalcFlux(true);

      OutputFlux(flux_);
    }

    // write mean values of scalar(s)
    OutputMeanScalars();

    // output of electrode status to screen and file (only if existing)
    OutputElectrodeInfo();

    // magnetic field (if existing)
    if (magneticfield_ != Teuchos::null)
      output_->WriteVector("magnetic_field", magneticfield_,IO::DiscretizationWriter::nodevector);
  }

  // NOTE:
  // statistics output for normal fluxes at boundaries was already done during Update()

  return;
} // ScaTraTimIntImpl::Output


/*==========================================================================*/
// scalar degrees of freedom and related
/*==========================================================================*/

/*----------------------------------------------------------------------*
 |  set initial field for phi                                 gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetInitialField(
    const INPAR::SCATRA::InitialField init,
    const int startfuncno)
{
  switch(init)
  {
  case INPAR::SCATRA::initfield_zero_field:
  {
    phin_-> PutScalar(0.0);
    phinp_-> PutScalar(0.0);
    break;
  }
  case INPAR::SCATRA::initfield_field_by_function:
  case INPAR::SCATRA::initfield_disturbed_field_by_function:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(0,lnode);

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function
        double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(k,lnode->X(),time_,NULL);
        int err = phin_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }

    // for NURBS discretizations we have to solve a least squares problem,
    // with high accuracy! (do nothing for Lagrangian polynomials)
    const Teuchos::ParameterList& scatradyn =
      DRT::Problem::Instance()->ScalarTransportDynamicParams();
    const int lstsolver = scatradyn.get<int>("LINEAR_SOLVER");
    if (lstsolver == (-1))
      dserror("no linear solver defined for least square NURBS problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number! Note: this solver block is misused for the least square problem. Maybe one should add a separate parameter for this.");

    DRT::NURBS::apply_nurbs_initial_condition(
        *discret_  ,
        errfile_,
        DRT::Problem::Instance()->SolverParams(lstsolver),
        startfuncno,
        phin_     );

    // initialize also the solution vector. These values are a pretty good guess for the
    // solution after the first time step (much better than starting with a zero vector)
    phinp_->Update(1.0,*phin_ ,0.0);

    // add random perturbation for initial field of turbulent flows
    if(init==INPAR::SCATRA::initfield_disturbed_field_by_function)
    {
      int err = 0;

      // random noise is relative to difference of max-min values of initial profile
      double perc = extraparams_->sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST",0.1);

      // out to screen
      if (myrank_==0)
      {
        cout << "Disturbed initial scalar profile:   max. " << perc*100 << "% random perturbation\n";
        cout << "\n\n";
      }

      // get overall max and min values and range between min and max
      double maxphi(0.0);
      double minphi(0.0);
      err = phinp_->MaxValue(&maxphi);
      if (err > 0) dserror("Error during evaluation of maximum value.");
      err = phinp_->MinValue(&minphi);
      if (err > 0) dserror("Error during evaluation of minimum value.");
      double range = abs(maxphi - minphi);

      // disturb initial field for all degrees of freedom
      for (int k=0; k < phinp_->MyLength(); ++k)
      {
        double randomnumber = DRT::Problem::Instance()->Random()->Uni();
        double noise = perc * range * randomnumber;
        err += phinp_->SumIntoMyValues(1,&noise,&k);
        err += phin_ ->SumIntoMyValues(1,&noise,&k);
        if (err!=0) dserror("Error while disturbing initial field.");
      }
    }
    break;
  }
  case INPAR::SCATRA::initfield_field_by_condition:
  {
    // set initial field for ALL existing scatra fields
    const string field = "ScaTra";
    const int numdofpernode = discret_->NumDof(discret_->lRowNode(0));
    vector<int> localdofs(numdofpernode);

    for (int i = 0; i < numdofpernode; i++)
    {
      localdofs[i] = i;
    }
    discret_->EvaluateInitialField(field,phin_,localdofs);

    // initialize also the solution vector. These values are a pretty good guess for the
    // solution after the first time step (much better than starting with a zero vector)
    phinp_->Update(1.0,*phin_ ,0.0);

    break;
  }
  // discontinuous 0-1 field for progress variable in 1-D
  case INPAR::SCATRA::initfield_discontprogvar_1D:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get coordinate
      const double x = lnode->X()[0];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        double initialval = 0.0;
        if (x > -EPS10) initialval = 1.0;

        int err = 0;
        err += phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
    break;
  }
  // reconstructed initial profile for progress variable in x2-direction from
  // Lessani and Papalexandris (2006), also used in Moureau et al. (2007, 2009),
  // for two-dimensional flame-vortex interaction problem (x2=0-200)
  case INPAR::SCATRA::initfield_flame_vortex_interaction:
  {
    // locations separating region 1 from region 2 and region 2 from region 3
    const double loc12 = 98.5;
    const double loc23 = 103.0;

    // define parameters for region 1 (exponential function for curve fitting)
    const double beta1  = 1.65;
    const double delta1 = 1.0;
    const double trans1 = 100.0;

    // define parameters for region 2 (linear function for curve fitting)
    const double abs2 = 0.0879;
    const double fac2 = 0.139309333;
    const double trans2 = 98.5;

    // define parameters for region 3 (exponential function for curve fitting)
    const double beta3  = 3.506209;
    const double delta3 = 4.28875;
    const double trans3 = 103.0;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // define variable
    double initialval = 0.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x2-coordinate
      const double x2 = lnode->X()[1];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        if (x2 < loc12-EPS10)
          initialval = (1.0-(1.0/beta1))*exp((x2-trans1)/delta1);
        else if (x2 > loc23+EPS10)
          initialval = 1.0-(exp((1.0-beta3)*(x2-trans3)/delta3)/beta3);
        else
          initialval = fac2*(x2-trans2) + abs2;

        int err = 0;
        err += phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
    break;
  }
  // initial mixture-fraction profile for Rayleigh-Taylor instability
  case INPAR::SCATRA::initfield_raytaymixfrac:
  {
    // define interface thickness, sinusoidal disturbance wave amplitude and pi
    const double delta = 0.002;
    const double alpha = 0.001;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x1- and x2-coordinate
      const double x1 = lnode->X()[0];
      const double x2 = lnode->X()[1];

      // interface disturbance
      //double x2_int = 0.05*cos(pi*(x1+0.5));
      //double x2_int = 0.05*cos(2.0*pi*x1);
      double x2_int = 0.0;
      x2_int -= cos(4*M_PI*x1);
      x2_int -= cos(14*M_PI*x1);
      x2_int -= cos(23*M_PI*x1);
      x2_int -= cos(28*M_PI*x1);
      x2_int -= cos(33*M_PI*x1);
      x2_int -= cos(42*M_PI*x1);
      x2_int -= cos(51*M_PI*x1);
      x2_int -= cos(59*M_PI*x1);
      x2_int *= alpha;

      const double value = (x2_int-x2)/(2.0*delta);

      // values required for tanh-distribution
      const double vp = exp(value);
      const double vm = exp(-value);

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        // compute tanh-distribution
        double initialval = 0.0;
        initialval = 0.5*(1.0+(vp-vm)/(vp+vm));

        int err = 0;
        err += phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
    break;
  }
  // initial field for skew convection of L-shaped domain
  case INPAR::SCATRA::initfield_Lshapeddomain:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x1- and x2-coordinate
      const double x1 = lnode->X()[0];
      const double x2 = lnode->X()[1];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);

        // compute initial values 0.0 or 1.0 depending on geometrical location
        double initialval = 0.0;
        if ((x1 <= 0.25 and x2 <= 0.5) or (x1 <= 0.5 and x2 <= 0.25))
          initialval = 1.0;

        int err = 0;
        err += phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good
        // guess for the solution after the first time step (much better than
        // starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
    break;
  }
  case INPAR::SCATRA::initfield_facing_flame_fronts:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x1- and x2-coordinate
      const double x1 = lnode->X()[0];
      //const double x2 = lnode->X()[1];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function

        double initialval;
        if (x1 < 0.0) initialval = -(x1+0.75);
        else initialval = x1-0.75;

        int err = 0;
        err += phin_->ReplaceMyValues(1,&initialval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
    break;
  }
  case INPAR::SCATRA::initfield_oracles_flame:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    const double eps = 0.00152;
    //const double xsing = 0.2;
    //const double zsing = 0.7525-0.05;//0.0354;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // get x1, x2 and x3-coordinate
      //const double x1 = lnode->X()[0];
      const double x2 = lnode->X()[1];
      //const double x3 = lnode->X()[2];

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function

        double initval = 0.0;

        // initial plane implementation for periodic spanwise boundary
        if (x2 >= 0.0)
          initval = (x2-0.0354) - eps;
        else
          initval = (-0.0354-x2) - eps;

#if 0
        // initial wedge implementation for periodic spanwise boundary
        if (x1 <= 0.0)
        {
          if (x2 >= 0.0)
            initval = (x2-0.0354) - eps;
          else
            initval = (-0.0354-x2) - eps;
        }
        else if (x1 > 0.0 and x1 < xsing)
        {
          initval = abs(x2)-0.0354*(xsing-x1)/xsing - eps;
        }
        else if (x1 >= xsing)
          initval = x1 - xsing - eps;
        else
          dserror("impossible!");
#endif

#if 0
        // initial wedge implementation for spanwise walls
        if (x1 <= 0.0)
        {
          if ( x3 <= -zsing and abs(x2) <= abs(x3+zsing) )
          {
            initval = (-0.7525-x3) - eps;
          }
          else if ( x3 >= zsing and abs(x2) <= (x3-zsing) )
          {
            initval = (x3-0.7525) - eps;
          }
          else if ( x2 >= 0.0 and ( x2 > abs(x3+zsing) or x2 > (x3-zsing) ))
          {
            initval = (x2-0.0354) - eps;
          }
          else if ( x2 < 0.0 and (-x2 > abs(x3+zsing) or -x2 > (x3-zsing) ))
          {
            initval = (-0.0354-x2) - eps;
          }
          else
            dserror("coordinate out of range of ORACLES initial function");
        }
        else if (x1 > 0.0 and x1 < xsing)
        {
          if (abs(x3) <= 0.07)
            initval = abs(x2)-0.0354*(xsing-x1)/xsing - eps;
          else
          {
            initval = 0.07525-0.07;
          }
        }
        else if (x1 >= xsing)
          initval = x1 - xsing - eps;
        else
          dserror("impossible!");
#endif
        int err = 0;
        err += phin_->ReplaceMyValues(1,&initval,&doflid);
        // initialize also the solution vector. These values are a pretty good guess for the
        // solution after the first time step (much better than starting with a zero vector)
        err += phinp_->ReplaceMyValues(1,&initval,&doflid);
        if (err != 0) dserror("dof not on proc");
      }
    }
  break;
  }
  default:
    dserror("Unknown option for initial field: %d", init); break;
  } // switch(init)

  return;
} // ScaTraTimIntImpl::SetInitialField


/*----------------------------------------------------------------------*
 | iterative update of concentrations                                   |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
{
  // store incremental vector to be available for convergence check
  // if incremental vector is received from outside for coupled problem
  increment_->Update(1.0,*inc,0.0);

  // update scalar values by adding increments
  phinp_->Update(1.0,*inc,1.0);
} // UpdateIter

/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
// general framework
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | prepare Krylov space projection                            gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareKrylovSpaceProjection()
{
  if (project_)
  {
    vector<DRT::Condition*> KSPcond;
    discret_->GetCondition("KrylovSpaceProjection",KSPcond);
    int nummodes = KSPcond.size();

    bool justcreated(false);
    // create vectors if not existing yet
    if (w_ == Teuchos::null)
    {
      w_ = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofRowMap()),nummodes,true));
      justcreated = true;
    }
    if (c_ == Teuchos::null)
    {
      c_ = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofRowMap()),nummodes,true));
      justcreated = true;
    }

    if (isale_ or justcreated) // fixed grid: compute w_,c_ only once at beginning!
    {
      for (int imode = 0; imode < nummodes; ++imode)
      {
        // zero w and c completely
        if (imode == 0)
        {
          w_->PutScalar(0.0);
          c_->PutScalar(0.0);
        }

        // in this case, we want to project out some zero pressure modes
        const string* definition = KSPcond[imode]->Get<string>("weight vector definition");

        // get rigid body modes
        const std::vector<double>* mode = KSPcond[imode]->Get<std::vector<double> >("mode");

        int numdof = 0;
        Epetra_IntSerialDenseVector dofids(6);
        for(int rr=0;rr<6;rr++)
        {
          if(abs((*mode)[rr])>1e-14)
          {
            numdof++;
            dofids(rr)=rr;
          }
          else
            dofids(rr)=-1;
        }

        if(*definition == "pointvalues")
        {
          dserror("option pointvalues not implemented");
        }
        else if(*definition == "integration")
        {
          Teuchos::ParameterList mode_params;

          // set parameters for elements
          mode_params.set<int>("action",SCATRA::integrate_shape_functions);
          mode_params.set<int>("scatratype",scatratype_);
          mode_params.set("dofids",dofids);

          mode_params.set("isale",isale_);
          if (isale_)
            AddMultiVectorToParameterList(mode_params,"dispnp",dispnp_);

          /* evaluate KrylovSpaceProjection condition in order to get
    // integrated nodal basis functions w_
    // Note that in the case of definition integration based,
    // the average pressure will vanish in an integral sense
    //
    //                    /              /                      /
    //   /    \          |              |  /          \        |  /    \
    //  | w_*p | = p_i * | N_i(x) dx =  | | N_i(x)*p_i | dx =  | | p(x) | dx = 0
    //   \    /          |              |  \          /        |  \    /
    //                   /              /                      /
           */

          // get an RCP of the current column Epetra_Vector of the MultiVector
          Teuchos::RCP<Epetra_Vector> wi = Teuchos::rcp((*w_)(imode),false);

          // compute integral of shape functions
          discret_->EvaluateCondition
              (mode_params           ,
              Teuchos::null      ,
              Teuchos::null      ,
              wi                 ,
              Teuchos::null      ,
              Teuchos::null      ,
              "KrylovSpaceProjection");

        }
        else
        {
          dserror("unknown definition of weight vector w for restriction of Krylov space");
        }

        // set the current kernel basis vector
        for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
        {
          DRT::Node* node = discret_->lRowNode(inode);
          vector<int> gdof = discret_->Dof(node);
          int numdof = gdof.size();
          if (numdof > 6) dserror("only up to 6 dof per node supported");
          for(int rr=0;rr<numdof;++rr)
          {
            const double val = (*mode)[rr];
            int err = c_->ReplaceGlobalValue(gdof[rr],imode,val);
            if (err != 0) dserror("error while inserting value into c_");
          }
        }

      } // loop over nummodes
    }
  }

  return;

} // ScaTraTimIntImpl::PrepareKrylovSpaceProjection

/*----------------------------------------------------------------------*
 | export multivector to column map & add it to parameter list gjb 06/09|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddMultiVectorToParameterList
(Teuchos::ParameterList& p,
    const std::string name,
    Teuchos::RCP<Epetra_MultiVector> vec
)
{
  if (vec != Teuchos::null)
  {
    //provide data in node-based multi-vector for usage on element level
    // -> export to column map is necessary for parallel evaluation
    //SetState cannot be used since this multi-vector is nodebased and not dofbased!
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    int numcol = vec->NumVectors();
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,numcol));
    LINALG::Export(*vec,*tmp);
    p.set(name,tmp);
  }
  else
    p.set(name,Teuchos::null);

  return;
} // SCATRA::ScaTraTimIntImpl::AddMultiVectorToParameterList

/*----------------------------------------------------------------------*
 | add approximation to flux vectors to a parameter list      gjb 05/10 |
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::AddFluxApproxToParameterList
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | application of Dirichlet boundary conditions                         |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletToSystem()
{
  // -------------------------------------------------------------------
  // Apply Dirichlet boundary conditions to system matrix
  // -------------------------------------------------------------------
  if (incremental_)
  {
    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the Dirichlet positions
    // are not used anyway.
    // We could avoid this though, if the dofrowmap would not include
    // the Dirichlet values as well. But it is expensive to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    //--------- Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    increment_->PutScalar(0.0);

    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }
  }
  else
  {
    // time measurement: application of DBC to system
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

    LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,phinp_,*(dbcmaps_->CondMap()));
  }
  return;
} // SCATRA::ScaTraTimIntImpl::ApplyDirichletToSystem

/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}           gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> phinp,
  Teuchos::RCP<Epetra_Vector> phidt
)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:      + apply dirich cond.");

  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time",time);  // actual time t_{n+1}

  // predicted Dirichlet values
  // \c  phinp then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->EvaluateDirichlet(p,phinp,phidt,Teuchos::null,Teuchos::null,dbcmaps_);
  discret_->ClearState();

  return;
} // SCATRA::ScaTraTimIntImpl::ApplyDirichletBC

/*----------------------------------------------------------------------*
 | compute outward pointing unit normal vectors at given b.c.  gjb 01/09|
 *----------------------------------------------------------------------*/
// RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::ComputeNormalVectors
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | evaluate Neumann inflow boundary condition                  vg 03/09 |
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::ComputeNeumannInflow
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | evaluate boundary cond. due to convective heat transfer     vg 10/11 |
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::EvaluateConvectiveHeatTransfer(
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | contains the residual scaling and addition of Neumann terms vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ScalingAndNeumann()
{
  // scaling to get true residual vector for all time integration schemes
  // in incremental case: boundary flux values can be computed from trueresidual
  if (incremental_) trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  AddNeumannToResidual();

  // add potential Neumann inflow or convective heat transfer boundary
  // conditions (simultaneous evaluation of both conditions not allowed!)
  if (neumanninflow_)     ComputeNeumannInflow(sysmat_,residual_);
  else if (convheatrans_) EvaluateConvectiveHeatTransfer(sysmat_,residual_);

  return;
} // ScaTraTimIntImpl::ScalingAndNeumann

/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions at t_{n+1}             gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyNeumannBC
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> phinp,
  Teuchos::RCP<Epetra_Vector> neumann_loads
)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // set time for evaluation of Neumann boundary conditions as parameter
  // depending on time-integration scheme
  Teuchos::ParameterList p;
  SetTimeForNeumannEvaluation(p);
  p.set<int>("scatratype",scatratype_);
  p.set("isale",isale_);
  // provide displacement field in case of ALE
  if (isale_) AddMultiVectorToParameterList(p,"dispnp",dispnp_);

  discret_->ClearState();
  // evaluate Neumann conditions at actual time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateNeumann(p,*neumann_loads);
  discret_->ClearState();

  return;
} // SCATRA::ScaTraTimIntImpl::ApplyNeumannBC

/*----------------------------------------------------------------------*
 | contains the assembly process for matrix and rhs            vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AssembleMatAndRHS()
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // zero out matrix entries
  sysmat_->Zero();

  // reset the residual vector
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  if (reinitswitch_ == true)
  {
    eleparams.set<int>("action",SCATRA::reinitialize_levelset);
  }
  else if(timealgo_ == INPAR::SCATRA::timeint_tg2
       or timealgo_ == INPAR::SCATRA::timeint_tg3)
  {
    // taylor galerkin transport of levelset
    eleparams.set<int>("action",SCATRA::calc_TG_mat_and_rhs);
  }
  else
  {
    // standard case
    eleparams.set<int>("action",SCATRA::calc_mat_and_rhs);
  }

  // DO THIS AT VERY FIRST!!!
  // compute reconstructed diffusive fluxes for better consistency
  const enum INPAR::SCATRA::Consistency consistency
  = DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(params_->sublist("STABILIZATION"),"CONSISTENCY");
  if (consistency == INPAR::SCATRA::consistency_l2_projection_lumped)
  {
    // compute flux approximation and add it to the parameter list
    AddFluxApproxToParameterList(eleparams,INPAR::SCATRA::flux_diffusive_domain);
  }

  // set type of scalar transport problem
  eleparams.set<int>("scatratype",scatratype_);

  // other parameters that might be needed by the elements
  eleparams.set("time-step length",dta_);
  eleparams.set("incremental solver",incremental_);
  eleparams.set<int>("form of convective term",convform_);
  eleparams.set<int>("fs subgrid diffusivity",fssgd_);
  // set general parameters for turbulent flow
  // prepare dynamic Smagorinsky model if required,
  // i.e. calculate turbulent Prandtl number
  if ((timealgo_ == INPAR::SCATRA::timeint_gen_alpha or timealgo_ == INPAR::SCATRA::timeint_one_step_theta
      or timealgo_ == INPAR::SCATRA::timeint_bdf2) and reinitswitch_ == false)
    DynamicComputationOfCs();
  eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");
  // and set parameters for multifractal subgrid-scale modeling
  if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales)
    eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");
  eleparams.set("turbulent inflow",turbinflow_);
  eleparams.set("frt",frt_);// ELCH specific factor F/RT
  if (scatratype_ == INPAR::SCATRA::scatratype_loma)
    eleparams.set<bool>("update material",(&(extraparams_->sublist("LOMA")))->get<bool>("update material",false));

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
  AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);
  AddMultiVectorToParameterList(eleparams,"magnetic field",magneticfield_);
  // and provide fine-scale velocity for multifractal subgrid-scale modeling only
  if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales or fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
    AddMultiVectorToParameterList(eleparams,"fine-scale velocity field",fsvel_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set switch for reinitialization
  eleparams.set("reinitswitch",reinitswitch_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  // parameters for Elch/DiffCond formulation
  if(IsElch(scatratype_))
    eleparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

  // set vector values needed by elements
  discret_->ClearState();

  // AVM3 separation for incremental solver: get fine-scale part of scalar
  if (incremental_ and
      (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
   AVM3Separation();

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(eleparams);

  // add reinitialization specific time-integration parameters
  if (reinitswitch_) AddReinitializationParameters(eleparams);

  // call loop over elements (with or without subgrid-diffusivity(-scaling) vector)
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,subgrdiff_,Teuchos::null);
  else
    discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

//  (SystemMatrix()->EpetraMatrix())->Print(cout); // kn nis

  discret_->ClearState();

  //----------------------------------------------------------------------
  // apply weak Dirichlet boundary conditions
  //----------------------------------------------------------------------
  {
    Teuchos::ParameterList mhdbcparams;

    // set action for elements
    mhdbcparams.set<int>("action",SCATRA::bd_calc_weak_Dirichlet);
    mhdbcparams.set("incremental solver",incremental_);
    mhdbcparams.set("isale",isale_);

    mhdbcparams.set<int>("scatratype",INPAR::SCATRA::scatratype_condif);

    AddMultiVectorToParameterList(mhdbcparams,"convective velocity field",convel_);
    AddMultiVectorToParameterList(mhdbcparams,"velocity field",vel_);
    AddSpecificTimeIntegrationParameters(mhdbcparams);

    // evaluate all mixed hybrid Dirichlet boundary conditions
    discret_->EvaluateConditionUsingParentData
      (mhdbcparams          ,
       sysmat_              ,
       Teuchos::null        ,
       residual_            ,
       Teuchos::null        ,
       Teuchos::null        ,
       "LineWeakDirichlet");

    discret_->EvaluateConditionUsingParentData
      (mhdbcparams          ,
       sysmat_              ,
       Teuchos::null        ,
       residual_            ,
       Teuchos::null        ,
       Teuchos::null        ,
       "SurfaceWeakDirichlet");

    // clear state
    discret_->ClearState();
  }

  AssembleMatAndRHS_Boundary();


  // AVM3 scaling for non-incremental solver: scaling of normalized AVM3-based
  // fine-scale subgrid-diffusivity matrix by subgrid diffusivity
  if (not incremental_ and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    AVM3Scaling(eleparams);

  // finalize the complete matrix
  sysmat_->Complete();

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpuele;

  if (msht_!=INPAR::FLUID::no_meshtying)
  {
    meshtying_->PrepareMeshtyingSystem(sysmat_, residual_);
  }

  return;
} // ScaTraTimIntImpl::AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | contains the linear solver                                  vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::LinearSolve()
{
  // -------------------------------------------------------------------
  //                        output to screen
  // -------------------------------------------------------------------
  PrintTimeStepInfo();

  // -------------------------------------------------------------------
  //                     preparations for solve
  // -------------------------------------------------------------------
  PrepareLinearSolve();

  // -------------------------------------------------------------------
  // Solve system in incremental or non-incremental case
  // -------------------------------------------------------------------
  if (incremental_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

    // get cpu time
    const double tcpusolve=Teuchos::Time::wallTime();

    solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,true);

    // end time measurement for solver
    dtsolve_=Teuchos::Time::wallTime()-tcpusolve;

    //------------------------------------------------ update solution vector
    UpdateIter(increment_);

    //--------------------------------------------- compute norm of increment
    double incnorm_L2(0.0);
    double scalnorm_L2(0.0);
    increment_->Norm2(&incnorm_L2);
    phinp_    ->Norm2(&scalnorm_L2);

    if (myrank_ == 0)
    {
      printf("+-------------------------------+-------------+\n");
      {
        if (scalnorm_L2 > EPS10)
          printf("|  relative increment (L2 norm) | %10.3E  |",incnorm_L2/scalnorm_L2);
        else // prevent division by an almost zero value
          printf("|  absolute increment (L2 norm) | %10.3E  |\n",incnorm_L2);
      }
      printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
      printf("+-------------------------------+-------------+\n");
    }
  }
  else
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

    // get cpu time
    const double tcpusolve=Teuchos::Time::wallTime();

    solver_->Solve(sysmat_->EpetraOperator(),phinp_,residual_,true,true);

    // end time measurement for solver
    dtsolve_=Teuchos::Time::wallTime()-tcpusolve;

    if (myrank_==0)
      printf("Solvertype linear_full (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
  }

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  ComputeIntermediateValues();

  return;
} // ScaTraTimIntImpl::LinearSolve

/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                       gjb 09/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::NonlinearSolve()
{
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:   + nonlin. iteration/lin. solve");

  bool stopgalvanostat(false);
  gstatnumite_=1;
  while (!stopgalvanostat) // galvanostatic control (ELCH)
  {
  // out to screen
  if(reinitswitch_==false)
    PrintTimeStepInfo();
  else
    PrintPseudoTimeStepInfoReinit();

  // special preparations for multifractal subgrid-scale model
  if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales)
    RecomputeMeanCsgsB();

  if (myrank_ == 0)
  {
    IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+\n"
             << "|- step/max -|- tol      [norm] -|-- con-res ---|-- pot-res ---|-- con-inc ---|-- pot-inc ---|-- con-res-inf ---|" << IO::endl;
  }

  // ---------------------------------------------- nonlinear iteration
  //stop nonlinear iteration when both increment-norms are below this bound
  const double  ittol = params_->sublist("NONLINEAR").get<double>("CONVTOL");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = (DRT::INPUT::IntegralValue<int>(params_->sublist("NONLINEAR"),"ADAPTCONV"));
  const double adaptolbetter = params_->sublist("NONLINEAR").get<double>("ADAPTCONV_BETTER");
  const double abstolres = params_->sublist("NONLINEAR").get<double>("ABSTOLRES");
  double       actresidual(0.0);

  int   itnum = 0;
  int   itemax = params_->sublist("NONLINEAR").get<int>("ITEMAX");
  bool  stopnonliniter = false;

  // perform explicit predictor step (-> better starting point for nonlinear solver)
  const bool explpredictor = (DRT::INPUT::IntegralValue<int>(params_->sublist("NONLINEAR"),"EXPLPREDICT") == 1);
  if (explpredictor)
    ExplicitPredictor();

/*
  const int numdim = 3;
  //create output file name
  std::stringstream temp;
  temp<< DRT::Problem::Instance()->OutputControlFile()->FileName()<<".nonliniter_step"<<step_;
  string outname = temp.str();
  string probtype = DRT::Problem::Instance()->ProblemType();

  RCP<IO::OutputControl> myoutputcontrol = Teuchos::rcp(new IO::OutputControl(discret_->Comm(),probtype,"Polynomial","myinput",outname,numdim,0,1000));
  // create discretization writer with my own control settings
  RCP<IO::DiscretizationWriter> myoutput =
    Teuchos::rcp(new IO::DiscretizationWriter(discret_,myoutputcontrol));
  // write mesh at step 0
  myoutput->WriteMesh(0,0.0);
*/

  while (stopnonliniter==false)
  {

#ifdef VISUALIZE_ELEDATA_GMSH
    const bool screen_out = false;
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("SubgridVelocityScatra", 0, 5, screen_out, 0);
    std::ofstream gmshfilecontent(filename.c_str());//, ios_base::out | ios_base::app);
    gmshfilecontent << "View \" " << "SubgridVelocityScatra" << " \" {\n";
    gmshfilecontent.close();
#endif

    itnum++;

    // check for negative/zero concentration values (in case of ELCH only)
    CheckConcentrationValues(phinp_);

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and rhs and assemble
    // -------------------------------------------------------------------
    AssembleMatAndRHS();

    // -------------------------------------------------------------------
    // potential residual scaling and potential addition of Neumann terms
    // -------------------------------------------------------------------
    ScalingAndNeumann();

    // add contributions due to electrode kinetics conditions
    EvaluateElectrodeKinetics(sysmat_,residual_);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the Dirichlet positions
    // are not used anyway.
    // We could avoid this though, if the dofrowmap would not include
    // the Dirichlet values as well. But it is expensive to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    // abort nonlinear iteration if desired
    if (AbortNonlinIter(itnum,itemax,ittol,abstolres,actresidual))
       break;

    //--------- Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    increment_->PutScalar(0.0);

    // Apply Dirichlet boundary conditions to system matrix
    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }

    //------------------------------------------------solve
    {
      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // time measurement: call linear solver
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + call linear solver");

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        solver_->AdaptTolerance(ittol,actresidual,adaptolbetter);
      }

/*
      // matrix printing options (DEBUGGING!)
      RCP<LINALG::SparseMatrix> A = SystemMatrix();
      if (A != Teuchos::null)
      {
        // print to file in matlab format
        const std::string fname = "sparsematrix.mtl";
        LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
        // print to screen
        (A->EpetraMatrix())->Print(cout);
        // print sparsity pattern to file
        LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
      }
      else
      {
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = BlockSystemMatrix();
        const std::string fname = "sparsematrix.mtl";
        LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
      }
      */
      // ScaleLinearSystem();  // still experimental (gjb 04/10)

      PrepareKrylovSpaceProjection();

      if (msht_!=INPAR::FLUID::no_meshtying)
        meshtying_->SolveMeshtying(*solver_, sysmat_, increment_, residual_, itnum, w_, c_, project_);
      else
        solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,itnum==1,w_,c_,project_);

      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_=Teuchos::Time::wallTime()-tcpusolve;
    }

    //------------------------------------------------ update solution vector
    phinp_->Update(1.0,*increment_,1.0);

    //-------- update values at intermediate time steps (only for gen.-alpha)
    ComputeIntermediateValues();

    // iteration number (only after that data output is possible)
  /*
    myoutput->NewStep(itnum,itnum);
    myoutput->WriteVector("phinp", phinp_);
   */

  } // nonlinear iteration

  stopgalvanostat = ApplyGalvanostaticControl();
  } // galvanostatic control

#ifdef VISUALIZE_ELEDATA_GMSH
  const bool screen_out = false;
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("SubgridVelocityScatra", 0, 5, screen_out, 0);
  std::ofstream gmshfilecontent(filename.c_str(), ios_base::out | ios_base::app);
  gmshfilecontent << "};\n";
  gmshfilecontent.close();
#endif

  return;
} // ScaTraTimIntImpl::NonlinearSolve

/*----------------------------------------------------------------------*
 | check if to stop the nonlinear iteration                    gjb 09/08|
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::AbortNonlinIter(
    const int itnum,
    const int itemax,
    const double ittol,
    const double abstolres,
    double& actresidual)
{
  //----------------------------------------------------- compute norms
  double incconnorm_L2(0.0);
  double incpotnorm_L2(0.0);

  double connorm_L2(0.0);
  double potnorm_L2(0.0);

  double conresnorm(0.0);
  double potresnorm(0.0);

  double conresnorminf(0.0);

  if (IsElch(scatratype_))
  {
    Teuchos::RCP<Epetra_Vector> onlycon = splitter_->ExtractOtherVector(residual_);
    onlycon->Norm2(&conresnorm);

    splitter_->ExtractOtherVector(increment_,onlycon);
    onlycon->Norm2(&incconnorm_L2);

    splitter_->ExtractOtherVector(phinp_,onlycon);
    onlycon->Norm2(&connorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(residual_);
    onlypot->Norm2(&potresnorm);

    splitter_->ExtractCondVector(increment_,onlypot);
    onlypot->Norm2(&incpotnorm_L2);

    splitter_->ExtractCondVector(phinp_,onlypot);
    onlypot->Norm2(&potnorm_L2);
  }
  else
  {
    residual_ ->Norm2(&conresnorm);
    increment_->Norm2(&incconnorm_L2);
    phinp_    ->Norm2(&connorm_L2);
    residual_ ->NormInf(&conresnorminf);
  }

  // care for the case that nothing really happens in the concentration
  // or potential field
  if (connorm_L2 < 1e-5)
  {
    connorm_L2 = 1.0;
  }
  if (potnorm_L2 < 1e-5)
  {
    potnorm_L2 = 1.0;
  }

  // absolute tolerance for deciding if residual is (already) zero
  // prevents additional solver calls that will not improve the residual anymore

  //-------------------------------------------------- output to screen
  /* special case of very first iteration step:
      - solution increment is not yet available
      - do not perform a solver call when the initial residuals are < EPS14*/
  if (itnum == 1)
  {
    if (myrank_ == 0)
    {
      IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
               << std::setw(10) << std::setprecision(3) << std::scientific << conresnorm << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   |      --      |      --      | "
               << std::setw(10) << std::setprecision(3) << std::scientific << conresnorminf << "       | (      --     ,te="
               << std::setw(10) << std::setprecision(3) << std::scientific << dtele_ << ")" << IO::endl;
    }
    // abort iteration, when there's nothing more to do
    if ((conresnorm < abstolres) && (potresnorm < abstolres))
    {
      // print 'finish line'
      if (myrank_ == 0)
      {
        IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << IO::endl;
      }
      return true;
    }
  }
  /* ordinary case later iteration steps:
      - solution increment can be printed
      - convergence check should be done*/
  else
  {
    // print the screen info
    if (myrank_ == 0)
    {
      IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
               << std::setw(10) << std::setprecision(3) << std::scientific << conresnorm << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << incconnorm_L2/connorm_L2 << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << incpotnorm_L2/potnorm_L2 << "   | "
               << std::setw(10) << std::setprecision(3) << std::scientific << conresnorminf << "       | (ts="
               << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_ << ",te="
               << std::setw(10) << std::setprecision(3) << std::scientific << dtele_ << ")" << IO::endl;
    }

    // this is the convergence check
    // We always require at least one solve. We test the L_2-norm of the
    // current residual. Norm of residual is just printed for information
    if (conresnorm <= ittol and potresnorm <= ittol and
        incconnorm_L2/connorm_L2 <= ittol and incpotnorm_L2/potnorm_L2 <= ittol)
    {
      if (myrank_ == 0)
      {
        // print 'finish line'
        IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << IO::endl;
        // write info to error file
        if (errfile_!=NULL)
        {
          fprintf(errfile_,"elch solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
              itnum,itemax,ittol,conresnorm,potresnorm,
              incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
        }
      }
      // yes, we stop the iteration
      return true;
    }

    // abort iteration, when there's nothing more to do! -> more robustness
    if ((conresnorm < abstolres) && (potresnorm < abstolres))
    {
      // print 'finish line'
      if (myrank_ == 0)
      {
        IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl;
      }
      return true;
    }

    // if not yet converged go on...
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if ((itnum == itemax))
  {
    if (myrank_ == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");

      if (errfile_!=NULL)
      {
        fprintf(errfile_,"elch divergent solve:   %3d/%3d  tol=%10.3E[L_2 ]  cres=%10.3E  pres=%10.3E  cinc=%10.3E  pinc=%10.3E\n",
            itnum,itemax,ittol,conresnorm,potresnorm,
            incconnorm_L2/connorm_L2,incpotnorm_L2/potnorm_L2);
      }
    }
    // yes, we stop the iteration
    return true;
  }

  // return the maximum residual value -> used for adaptivity of linear solver tolerance
  actresidual = max(conresnorm,potresnorm);
  actresidual = max(actresidual,incconnorm_L2/connorm_L2);
  actresidual = max(actresidual,incpotnorm_L2/potnorm_L2);

  // check for INF's and NaN's before going on...
  if (std::isnan(incconnorm_L2) or
      std::isnan(incpotnorm_L2) or
      std::isnan(connorm_L2) or
      std::isnan(potnorm_L2) or
      std::isnan(conresnorm) or
      std::isnan(potresnorm))
    dserror("calculated vector norm is NaN.");

  if (abs(std::isinf(incconnorm_L2)) or
      abs(std::isinf(incpotnorm_L2))  or
      abs(std::isinf(connorm_L2))  or
      abs(std::isinf(potnorm_L2))  or
      abs(std::isinf(conresnorm))  or
      abs(std::isinf(potresnorm)) )
    dserror("calculated vector norm is INF.");

  return false;
} // ScaTraTimIntImpl::AbortNonlinIter

/*----------------------------------------------------------------------*
| returns matching string for each time integration scheme   gjb 08/08 |
*----------------------------------------------------------------------*/
std::string SCATRA::ScaTraTimIntImpl::MapTimIntEnumToString
(
   const enum INPAR::SCATRA::TimeIntegrationScheme term
)
{
  // length of return string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPAR::SCATRA::timeint_one_step_theta :
    return "One-Step-Theta";
    break;
  case INPAR::SCATRA::timeint_bdf2 :
    return "    BDF2      ";
    break;
  case INPAR::SCATRA::timeint_stationary :
    return "  Stationary  ";
    break;
  case INPAR::SCATRA::timeint_gen_alpha :
    return "  Gen. Alpha  ";
    break;
  case INPAR::SCATRA::timeint_tg2 :
    return "  Taylor Galerkin 2rd order  ";
    break;
  case INPAR::SCATRA::timeint_tg3 :
    return "  Taylor Galerkin 3rd order  ";
    break;
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }

  return "";
} // ScaTraTimIntImpl::MapTimIntEnumToString

/*----------------------------------------------------------------------*
 |  write current state to BINIO                             gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputState()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  // convective velocity (not written in case of coupled simulations)
//  if (cdvel_ != INPAR::SCATRA::velocity_Navier_Stokes)
//    output_->WriteVector("convec_velocity", convel_,IO::DiscretizationWriter::nodevector);

  // displacement field
  if (isale_) output_->WriteVector("dispnp", dispnp_);

  return;
} // ScaTraTimIntImpl::OutputState

/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        henke   12/09|
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::OutputToGmsh(
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::OutputFlux(RCP<Epetra_MultiVector> flux)
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | increment time and step for next iteration                     mr. x |
 *----------------------------------------------------------------------*/
inline void SCATRA::ScaTraTimIntImpl::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dta_;
}

/*----------------------------------------------------------------------*
 | time update of time-dependent materials                    gjb 07/12 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ElementMaterialTimeUpdate()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  p.set<int>("action", SCATRA::time_update_material);
  // further required parameters
  p.set<int>("scatratype",scatratype_);
  p.set("time-step length",dta_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // go to elements
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
      Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->ClearState();
  return;
}

/*==========================================================================*/
// ELCH
/*==========================================================================*/

// all defined in scalar_timint_implicit_service

/*==========================================================================*/
// AVM3
/*==========================================================================*/

// all defined in scalar_timint_implicit_service

/*==========================================================================*/
// functions used for reinitialization of level sets
/*==========================================================================*/

// all defined in scalar_timint_reinitialization

/*==========================================================================*/
//  obsolete or unused methods - to be deleted soon (at noon)!!!
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | scale lines of linear system prior to solve call           gjb 04/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ScaleLinearSystem()
{
  if (IsElch(scatratype_))
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = BlockSystemMatrix();
    if (A != Teuchos::null)
    {
    // scale ELCH Matrix:
    Teuchos::RCP<Epetra_Vector> d = splitter_->Vector(0);
    Teuchos::RCP<Epetra_Vector> scalefactors = splitter_->Vector(0);
    Teuchos::RCP<Epetra_Vector> scalefactorsENC = splitter_->Vector(1);
    scalefactorsENC->PutScalar(1.0);

    // access values located at main diagonal
    BlockSystemMatrix()->Matrix(0,0).ExtractDiagonalCopy(*d);
    double firstrowvalue(0.0);
    for (int r=0;r < d->MyLength(); r++)
    {
      double* values = d->Values();
      if (r % numscal_ == 0)
      {
        firstrowvalue = values[r];
        if (abs(firstrowvalue)<EPS14) dserror("diagonal value too small");
        int err = scalefactors->ReplaceMyValue(r,0,1.0); // set scalefactor 1.0
        if (err != 0) dserror("Error detected");
      }
      else
      {
        if(abs(values[r])<EPS14) dserror("devision by zero prevented");
        double scalefactor = 1.0*abs(firstrowvalue/values[r]); //sign-independent!
        int err = scalefactors->ReplaceMyValue(r,0,scalefactor);
        if (err != 0) dserror("Error detected");
      }
    } // for (int r=0;r < w->MyLength(); r++)
#if 0
    scalefactors->Print(cout);
    d->Print(cout);
    (((BlockSystemMatrix()->Matrix(0,0)).EpetraMatrix()->Print(cout)));
#endif
    // scale the complete rows and rhs
    int err = BlockSystemMatrix()->Matrix(0,0).LeftScale(*scalefactors);
       err += BlockSystemMatrix()->Matrix(0,1).LeftScale(*scalefactors);
// ENC scaling
       err += BlockSystemMatrix()->Matrix(1,0).LeftScale(*scalefactorsENC);

    Teuchos::RCP<Epetra_Vector> onlyconc = splitter_->ExtractOtherVector(residual_);
#if 0
    cout<<"residual:\n";
    residual_->Print(cout);
    cout<<"non-modified onlyconc:\n";
    onlyconc->Print(cout);
#endif
    //Multiply a Epetra_MultiVector with another, element-by-element.
    onlyconc->Multiply(1.0,*onlyconc,*scalefactors,0.0);
    // insert values into the whole rhs vector
    splitter_->InsertOtherVector(onlyconc,residual_);

    Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(residual_);
    //Multiply a Epetra_MultiVector with another, element-by-element.
    onlypot->Multiply(1.0,*onlypot,*scalefactorsENC,0.0);
    // insert values into the whole rhs vector
    splitter_->InsertCondVector(onlypot,residual_);

#if 0
    cout<<"modified onlyconc:\n";
    onlyconc->Print(cout);
    cout<<"modified residual:\n";
    residual_->Print(cout);
#endif
    if (err>0) dserror("Error during pre-scaling of linear system");


    {
      // matrix printing options (DEBUGGING!)
      RCP<LINALG::SparseMatrix> A = SystemMatrix();
      if (A != Teuchos::null)
      {
      }
      else
      {
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = BlockSystemMatrix();
        const std::string fname = "sparsematrix_scaled.mtl";
        LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
      }
    }

    if (myrank_==0)
    cout<<"Pre-Scaling of Linear System done."<<endl;

    //scalefactors->Print(cout);

    }
  } // pre-scale equation system for ELCH applications

  return;
} // ScaTraTimIntImpl::ScaleLinearSystem

/*----------------------------------------------------------------------*
 | construct toggle vector for Dirichlet dofs                  gjb 11/08|
 | assures backward compatibility for avm3 solver; should go away once  |
 *----------------------------------------------------------------------*/
// const Teuchos::RCP<const Epetra_Vector> SCATRA::ScaTraTimIntImpl::DirichletToggle()
// defined in scalar_timint_implicit_service.cpp

