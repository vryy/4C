/*----------------------------------------------------------------------------*/
/**
\file xcontact_algorithm_base.cpp

\brief base class of the inequality level-set approach algorithm for contact
       problems a.k.a. xcontact or extended contact

\maintainer Michael Hiermeier

\date Jun 15, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xcontact_algorithm_base.H"

#include "../drt_contact/contact_utils.H"

//#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_xfem/xfem_discretization_utils.H"

#include "../drt_io/io.H"

#include "../drt_inpar/inpar_scatra.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../drt_adapter/ad_str_xcontact.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_xstructure_algorithm.H"

#include "xcontact_multi_discretization_wrapper.H"
#include "xcontact_levelset_algorithm.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::ALGORITHM::Base::Base()
    : ADAPTER::AlgorithmBase(DRT::Problem::Instance()->GetDis("structure")->Comm(),
        DRT::Problem::Instance()->XContactDynamicParams()),
      isinit_(false),
      issetup_(false),
      p_xcontact_dyn_ptr_(Teuchos::null),
      p_struct_dyn_ptr_(Teuchos::null),
      p_scatra_dyn_ptr_(Teuchos::null),
      p_xfem_general_ptr_(Teuchos::null),
      structure_ptr_(Teuchos::null),
      scatra_ptr_(Teuchos::null),
      multi_discret_ptr_(Teuchos::null),
      i_scatra_contact_coupling_(Teuchos::null),
      numstep_(0),
      num_dof_per_node_(-1),
      max_num_dof_sets_(-1)
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::Init(
    const Teuchos::ParameterList& p_xcontact_dyn,
    const Teuchos::ParameterList& p_structure_dyn,
    const Teuchos::ParameterList& p_scatra_dyn,
    const Teuchos::ParameterList& p_xfem_general,
    const Teuchos::RCP<DRT::Discretization>& structdis_ptr)
{
  issetup_ = false;

  // copy and store the XCONTACT DYNAMIC parameter list
  p_xcontact_dyn_ptr_ = Teuchos::rcp( new Teuchos::ParameterList( p_xcontact_dyn ) );
  p_struct_dyn_ptr_   = Teuchos::rcp( new Teuchos::ParameterList( p_structure_dyn ) );
  p_scatra_dyn_ptr_   = Teuchos::rcp( new Teuchos::ParameterList( p_scatra_dyn ) );
  p_xfem_general_ptr_ = Teuchos::rcp( new Teuchos::ParameterList( p_xfem_general ) );

  max_num_dof_sets_ = XFEMGeneralParams().get<int>("MAX_NUM_DOFSETS");

  // get the number of DoF's per node
  int gid_node = structdis_ptr->NodeRowMap()->MinMyGID();
  DRT::Node* node_ptr = structdis_ptr->gNode( gid_node );
  num_dof_per_node_ = structdis_ptr->NumDof( node_ptr );

  // check the fill complete status of the structure discretization
  if (not structdis_ptr->Filled())
    structdis_ptr->FillComplete();
  // store the structure discret pointer in the multi discretization wrapper
  multi_discret_ptr_ = Teuchos::rcp(new XCONTACT::MultiDiscretizationWrapper());
  multi_discret_ptr_->Init( "XContact_Wrapper", structdis_ptr,
        NumDofPerNode() * MaxNumDofSets() );

  /* specify remaining control parameters (which are not already set in the
   * ADAPTER::AlgorithmBase class.) */
  numstep_ = p_xcontact_dyn.get<int>( "NUMSTEP" );

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::Setup()
{
  CheckInit();

  // setup the structural side
  SetupStructure();

  // setup the scalar transport part
  SetupScaTra();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::SetupStructure()
{
  CheckInit();

  // get the contact condition groups
  std::vector<std::vector<DRT::Condition*> > ccond_grps(0);
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps,StructDiscret());
  if (ccond_grps.size()>1)
    dserror("Currently we support only ONE contact interface. This "
        "can be easily extended. But keep in mind, that each slave "
        "side needs a own ScaTra discretization!");

  /*------------------------------------------------------------------------*
   * split the structural discretization into a XFEM discretization, which
   * is located next to the (enriched) conditioned interfaces and a standard
   * discretization
   *
   *           ___boundary cond___    _____ enriched and conditioned boundary
   *          /                   \  /      interface nodes (o)
   *                                /
   *          o---o---o---o---o---o
   *         /   /   /   /   /   /|       enriched element row (xFem discret.)
   *        o---o---o---o---o---o +   <== (enriched (o) and std. nodes (+))
   *        | 0 | 1 | 2 | 3 | 4 |/|
   *        +---+---+---+---+---+ +       standard element row (std. discret.)
   *        | 5 | 6 | 7 | 8 | 9 |/    <== (only std. nodes (+))
   *        +---+---+---+---+---+
   *
   *                 __
   *                |  |
   *               _|  |_
   *               \    /
   *                \  /
   *                 \/
   *
   * We get a new cut xFem discretization, which is connected to the
   * conditioned boundary interface (o)
   *
   *          o---o---o---o---o---o
   *         /   /   /   /   /   /|
   *        o---o---o---o---o---o +
   *        | 0 | 1 | 2 | 3 | 4 |/   <== xstruct_dis_ptr
   *        +---+---+---+---+---+
   *
   * and the remaining standard structure discretization (+)
   *
   *          +---+---+---+---+---+
   *         /   /   /   /   /   /|
   *        +---+---+---+---+---+ +
   *        | 5 | 6 | 7 | 8 | 9 |/   <== struct_dis_ptr_
   *        +---+---+---+---+---+
   *
   * The two discretizations share the same node ID's at the coupling interface,
   * but differ in the global degrees of freedom ID's!
   *------------------------------------------------------------------------*/
  Teuchos::RCP<DRT::DiscretizationInterface> xstruct_dis_ptr =
      Teuchos::rcp(new DRT::DiscretizationXFEM("xstructure",
          Teuchos::rcp(StructDiscret().Comm().Clone())));
  xstruct_dis_ptr->FillComplete(false,false,false);
  const XFEM::UTILS::XFEMDiscretizationBuilder xdis_builder;
  xdis_builder.SetupXFEMDiscretization(XFEMGeneralParams(),
      StructDiscretPtr(),xstruct_dis_ptr,ccond_grps[0]);

  // --------------------------------------------------------------------------
  // add the xfem discretization to the multi discretization wrapper
  // --------------------------------------------------------------------------
  multi_discret_ptr_->AddDiscretization(xstruct_dis_ptr);
  multi_discret_ptr_->Setup();

  // --------------------------------------------------------------------------
  // show wrapped discretizations
  // --------------------------------------------------------------------------
  multi_discret_ptr_->Print(std::cout);

  // --------------------------------------------------------------------------
  // Build structure algorithm
  // --------------------------------------------------------------------------
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> structure_alg_base =
      Teuchos::rcp(new ADAPTER::XStructureAlgorithm());

  // call init and setup
  structure_alg_base->Init(XContactDynParams(),StructDynParams(),
      multi_discret_ptr_);
  structure_alg_base->Setup();

  // get the actual field pointer
  structure_ptr_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureXContact>(
      structure_alg_base->StructureField(),true);
  structure_ptr_->Init(XFEMGeneralParams(),NumDofPerNode());
  structure_ptr_->Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::SetupScaTra()
{
  CheckInit();

  // ------------------------------------------------------------------------
  // create new ScaTra discretizations
  // ------------------------------------------------------------------------
  XCONTACT::MultiDiscretizationWrapper::XDisPairedPtrVector::const_iterator cit;
  for (cit  = multi_discret_ptr_->GetContactIDiscret().begin();
       cit != multi_discret_ptr_->GetContactIDiscret().end(); ++cit)
  {
    Teuchos::RCP<DRT::Discretization> scatra_dis_ptr =
        Teuchos::rcp(new DRT::Discretization("scatra",
            Teuchos::rcp<Epetra_Comm>(StructDiscret().Comm().Clone())));
    scatra_dis_ptr->FillComplete(false,false,false);

    DRT::DiscretizationInterface & discret =
        multi_discret_ptr_->Discret(cit->first);
    std::vector<DRT::Condition*> conds(1,ExtractSlaveCondition(discret));

    // clone the wrapped contact discretization
    DRT::UTILS::CloneDiscretizationFromCondition<SCATRA::ScatraFluidCloneStrategy>(
        discret,*scatra_dis_ptr,conds);

    // give scatra new dof-sets (starts after structure)
    Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::DofSet());
    scatra_dis_ptr->ReplaceDofSet(newdofset,true);
    scatra_dis_ptr->FillComplete(true,true,true);

    // set the level-set implementation type
    for(int i=0; i<scatra_dis_ptr->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(
          scatra_dis_ptr->lColElement(i));

      // set the implementation type of the level-set elements
      if(element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(INPAR::SCATRA::impltype_lsreinit);
    }

    // create discretization writer
    scatra_dis_ptr->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatra_dis_ptr)));

    scatra_dis_ptr->Print(std::cout);

    DRT::Problem::Instance()->AddDis("scatra",scatra_dis_ptr);
    multi_discret_ptr_->AddScaTraIDiscret(cit->first, scatra_dis_ptr);
  }

  /* Setup the contact scatra dof exporters here once (in the default case
   * this is done during the FillComplete call, but we don't want to
   * rebuild everything once again. */
  multi_discret_ptr_->SetupContactScaTraDofExporter();

  // --------------------------------------------------------------------------
  // Build scalar transport algorithm
  // --------------------------------------------------------------------------
  SetInitialScaTraParams(ScaTraDynParams());

  // get linear solver id from SCALAR TRANSPORT DYNAMIC
  const int linsolvernumber = ScaTraDynParams().get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
   dserror("No linear solver defined for the eXtended contact problem. "
       "Please set LINEAR_SOLVER in the SCALAR TRANSPORT DYNAMIC section to a "
       "valid number!");
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_alg_base =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  // call init and setup
  scatra_alg_base->Init(XContactDynParams(),
          ScaTraDynParams(),DRT::Problem::Instance()->SolverParams(linsolvernumber),
          "scatra",false);
  scatra_alg_base->Setup();

  // get the actual field pointer
  scatra_ptr_ = Teuchos::rcp_dynamic_cast<XCONTACT::LEVELSET::Algorithm>(
      scatra_alg_base->ScaTraField(),true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Condition * XCONTACT::ALGORITHM::Base::ExtractSlaveCondition(
    const DRT::DiscretizationInterface & discret) const
{
  std::vector<std::vector<DRT::Condition*> > ccond_grps(0);
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps,discret);

  if (ccond_grps.size() > 1)
    dserror("Something went wrong. The wrapped discretization %s is supposed "
        "to hold only one contact condition group!",
        discret.Name().c_str());

  // --- Extract the slave condition from the structural discretization -------
  std::vector<bool> isslave(0);
  std::vector<bool> isself(0);
  CONTACT::UTILS::GetMasterSlaveSideInfo(isslave,isself,ccond_grps[0]);

  if (isself[0])
    dserror("Self contact is currently unsupported!");

  // get the condition from the slave side
  DRT::Condition* slave_cond = NULL;
  for (std::size_t i=0;i<isslave.size();++i)
  {
    if (isslave[i])
    {
      slave_cond = ccond_grps[0][i];
      break;
    }
  }

  if (slave_cond==NULL)
    dserror("The slave side condition could not be extracted!");

  return slave_cond;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::Timeloop()
{
  OutputInitialField();

  // time loop
  while (NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop for particular type of algorithm
    OuterLoop();

    // prepare output
    PrepareOutput();

    // update for next time step
    Update();

    // write output to screen and files
    Output();

  } // time loop

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::OutputInitialField()
{
  // initial output of the structural field
  StructureField().PrepareOutput();
  StructureField().Output();

  /* initial output of the scaTra field (should do nothing at the moment,
   * since an initial contact situation is not supported!) */
  ScaTraField().Output();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::PrepareTimeStep()
{
  IncrementTimeAndStep();

  // prepare structure time step, among other things, predict displacement field
  StructureField().PrepareTimeStep();

  // prepare scalar transport time step
  // (+ initial signed distance function will be evaluated, if the bodies
  //    come into contact)
  ScaTraField().PrepareTimeStep();

  // Synchronism check between structure and level-set algorithm
  if (StructureField().Time() != Time())
    dserror("Time in Structure time integration differs from time in "
        "eXtended contact algorithm \n(t_structure = %d, t_xcontact_algo = %d",
        StructureField().Time(),Time());

  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra time integration differs from time in "
        "eXtended contact algorithm \n(t_scatra = %d, t_xcontact_algo = %d",
        ScaTraField().Time(),Time());;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::Update()
{
  // update the scalar transport field
  ScaTraField().Update();

  // update the structure field and the time step!
  StructureField().Update();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::PrepareOutput()
{
  StructureField().PrepareOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::Output()
{
  /* NOTE: The order is important here! Herein, control file entries are
   * written, defining the order in which the filters handle the
   * discretizations, which in turn defines the dof number ordering of the
   * discretizations. */
  StructureField().Output();
  ScaTraField().Output();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::SetScaTraValuesInStructure()
{
  StructureField().SetScaTraValuesInStructure_Np( *ScaTraField().Phinp() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XCONTACT::ALGORITHM::Base::SetStructureValuesInScaTra()
{
  // get the contact interface velocity
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<DRT::DiscretizationInterface>& XCONTACT::ALGORITHM::Base::
    StructDiscretPtr()
{
  return multi_discret_ptr_->DiscretPtr( XFEM::structure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::DiscretizationInterface& XCONTACT::ALGORITHM::Base::StructDiscret()
{
  CheckInit();
  return multi_discret_ptr_->Discret( XFEM::structure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<DRT::DiscretizationInterface>& XCONTACT::ALGORITHM::Base::
    XStructDiscretPtr()
{
  return multi_discret_ptr_->DiscretPtr( XFEM::xstructure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::DiscretizationInterface& XCONTACT::ALGORITHM::Base::XStructDiscret()
{
  CheckInit();
  return multi_discret_ptr_->Discret( XFEM::xstructure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<DRT::Discretization>& XCONTACT::ALGORITHM::Base::ScaTraDiscretPtr()
{
  CheckInit();
  return multi_discret_ptr_->ScaTraDiscretPtr( XFEM::xstructure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& XCONTACT::ALGORITHM::Base::ScaTraDiscret()
{
  CheckInit();
  return multi_discret_ptr_->ScaTraDiscret( XFEM::xstructure );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Base::SetInitialScaTraParams(
    Teuchos::ParameterList& p_scatra_dyn) const
{
  // skip the initial time derivative
  p_scatra_dyn.set<std::string>("SKIPINITDER", "Yes");
}
