/*!----------------------------------------------------------------------

\brief specialization of ssi2wc, including "structale"-surface growth

\level 2

\maintainer Jonas Eichinger
 *----------------------------------------------------------------------*/
#include "ssi_partitioned_2wc_protrusionformation.H"
#include "../drt_inpar/inpar_cell.H"

#include "../drt_adapter/ad_ale_fluid.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structalewrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_scatra_wrapper_cellmigration.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex20.H"
#include "../drt_so3/so_hex27.H"
#include "../drt_so3/so_tet4.H"
#include "../drt_so3/so_tet10.H"
#include "../drt_so3/so3_scatra.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_base.H"

#include "../drt_io/io_control.H"
#include "../drt_wear/wear_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_clonestrategy.H"


/*----------------------------------------------------------------------*
 | constructor                                              rauch 01/16 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_PROTRUSIONFORMATION::SSI_Part2WC_PROTRUSIONFORMATION(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC(comm, globaltimeparams),
      omega_(-1.0),
      myrank_(comm.MyPID()),
      numdof_actin_(-1),
      ale_fluid_wrapper_(Teuchos::null),
      ssi_wrapper_(Teuchos::null),
      struct_ale_wrapper_(Teuchos::null),
      cell_scatra_wrapper_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------*
 | Initialize this object                                   rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_PROTRUSIONFORMATION::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part2WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // check if scatra in cell is set up with ale description
  if (not ScaTraField()->ScaTraField()->IsALE())
    dserror("We need an ALE description for the cell-scatra field!");

  return returnvar;
}


/*----------------------------------------------------------------------*
 | Init this object                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::Setup()
{
  // call init in base class
  SSI::SSI_Part2WC::Setup();

  // safety check
  if (StructAle() == Teuchos::null) dserror("No structure wrapper for struct-ale set!");

  // initialize pointer to structure
  ssi_wrapper_ = Teuchos::rcp_dynamic_cast<ADAPTER::SSIStructureWrapper>(StructureField());
  if (ssi_wrapper_ == Teuchos::null)
    dserror(
        "cast from ADAPTER::Structure to ADAPTER::MultiphysicsStructureWrapperCellMigration "
        "failed");

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(new ADAPTER::AleBaseAlgorithm(
      DRT::Problem::Instance()->SSIControlParams(), DRT::Problem::Instance()->GetDis("ale")));
  ale_fluid_wrapper_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFluidWrapper>(ale->AleField());
  if (ale_fluid_wrapper_ == Teuchos::null)
    dserror("cast from ADAPTER::Ale to ADAPTER::AleFsiWrapper failed");
  // create empty operator
  AleField()->CreateSystemMatrix();

  // build coupling objects for dof transfer between structure and ale
  const int ndim = DRT::Problem::Instance()->NDim();

  // create and seup coupling objects for matching grid
  coupalestru_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupalestru_->SetupCoupling(*AleField()->Discretization(), *StructureField()->Discretization(),
      *(StructureField()->Discretization()->NodeRowMap()),
      *(AleField()->Discretization()->NodeRowMap()), ndim, true, 1e-06);

  // set-up aux FSI interface
  Teuchos::RCP<STR::AUX::MapExtractor> interface = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface->Setup(*(scatra_->ScaTraField()->Discretization()),
      *(scatra_->ScaTraField()->Discretization()->DofRowMap()));

  // growth state
  x_ = Teuchos::rcp(
      new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth from last time step
  growth_n_ = Teuchos::rcp(
      new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth at new time step
  growth_np_ = Teuchos::rcp(
      new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth step increment growth_{n+1} = growth_{n} + growth_{step}
  growth_step_ = Teuchos::rcp(
      new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // initialize growth increment vector
  growthinc_ = Teuchos::rcp(
      new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true),

  // initialize delta_ale_ vector
      delta_ale_ = Teuchos::rcp(new Epetra_Vector(AleField()->Dispnp()->Map(), true));

  // source/sink vector to be added to scatra rhs
  sources_ = LINALG::CreateVector(*(scatra_->ScaTraField()->Discretization()->DofRowMap()), true);

  // numdof actin
  numdof_actin_ = DRT::Problem::Instance()
                      ->CellMigrationParams()
                      .sublist("PROTRUSION MODULE")
                      .get<int>("NUMDOF_ACTIN");
  int numdof_ARP23 = DRT::Problem::Instance()
                         ->CellMigrationParams()
                         .sublist("PROTRUSION MODULE")
                         .get<int>("NUMDOF_BRANCHES");
  int numdof_BarbedEnds = DRT::Problem::Instance()
                              ->CellMigrationParams()
                              .sublist("PROTRUSION MODULE")
                              .get<int>("NUMDOF_BARBEDENDS");
  if (myrank_ == 0)
  {
    std::cout << "\n  Number of Actin Monomer Dof (Volume) in Scalar Transport System: "
              << numdof_actin_ << std::endl;
    std::cout << "  Number of ARP2/3 Dof (Volume) in Scalar Transport System: " << numdof_ARP23
              << std::endl;
    std::cout << "  Number of Barbed End Dof (Volume) in Scalar Transport System: "
              << numdof_BarbedEnds << "\n"
              << std::endl;
  }

  // relaxation
  omega_ = DRT::Problem::Instance()
               ->CellMigrationParams()
               .sublist("PROTRUSION MODULE")
               .get<double>("RELAX_GROWTH");
  if (myrank_ == 0)
    std::cout << "\n  Use fixed relaxation parameter for protrusion formation omega = " << omega_
              << "\n"
              << std::endl;

  // build map todo replace FSICoupling by cell specific condition
  BuildConditionDofRowMap((StructureField()->Discretization())->GetCondition("FSICoupling"),
      StructureField()->Discretization(), conditiondofrowmap_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_PROTRUSIONFORMATION::InitFieldCoupling(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& scatra_disname)
{
  int returnvar = 0;

  // call SetupDiscretizations in base class
  returnvar = SSI::SSI_Base::InitFieldCoupling(comm, struct_disname, scatra_disname);

  return returnvar;
}


/*----------------------------------------------------------------------*
 |                                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::InitDiscretizations(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& scatra_disname)
{
  // call setup in base class
  SSI::SSI_Part2WC::InitDiscretizations(comm, struct_disname, scatra_disname);

  // new ale part
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);

  // clone ale from structure for ssi-actin assembly
  // access the ale discretization
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = DRT::Problem::Instance()->GetDis("ale");
  if (!aledis->Filled()) aledis->FillComplete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->NumGlobalNodes() == 0) dserror("ERROR: Structure discretization is empty!");

  // clone ale mesh from structure discretization
  if (aledis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(structdis, aledis);
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else
    dserror(
        "ERROR: Reading an ALE mesh from the input file is not supported for this problem type.\n"
        "NumGlobalNodes=%d",
        aledis->NumGlobalNodes());

  return;
}


/*----------------------------------------------------------------------*
 | Solve structure filed                                    rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  do structure step
  // -------------------------------------------------------------------
  StructureField()->Solve();

  // -------------------------------------------------------------------
  //                  do growth (ale) step
  // -------------------------------------------------------------------
  DoAleStep(x_);

  // application of mesh displacements to structural field,
  // update material displacements
  UpdateMatConf();

  // update dispnp and velnp
  UpdateSpatConf();

  // set mesh displacement and velocity fields
  SetStructSolution(structure_->Dispnp(), structure_->Velnp());

  // -------------------------------------------------------------------
  //                  evaluate new sources
  // -------------------------------------------------------------------
  return EvaluateSources();
}


/*----------------------------------------------------------------------*
 | Solve ale field                                          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::DoAleStep(Teuchos::RCP<Epetra_Vector> growthincrement)
{
  // initialize norm
  double normofgrowth = -1234.0;
  // calc norm
  growthincrement->Norm2(&normofgrowth);

  if (myrank_ == 0)
    std::cout << "=======================================\n"
                 "DoAleStep\n"
                 "Norm of applied growth "
              << std::setprecision(10) << normofgrowth << "\n"
              << "=======================================" << std::endl;

  // get lagrangian structure displacements
  Teuchos::RCP<Epetra_Vector> dispnp_ale = StructureToAle(StructureField()->Dispnp());

  // update ale field with lagrangian structure displacements
  AleField()->WriteAccessDispnp()->Update(1.0, *dispnp_ale, 0.0);

  // application of interface displacements as dirichlet conditions
  AleField()->ApplyInterfaceDisplacements(growthincrement);

  // solve time step todo remove auxiliary use of fsi stuff
  AleField()->TimeStep(ALE::UTILS::MapExtractor::dbc_set_part_fsi);

  return;
}


/*----------------------------------------------------------------------*
 | Solve Scatra field                                       rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  if (CellScatra() == Teuchos::null)
    dserror("Pointer to scatra wrapper for cell-scatra time integration is not set!");

  // rate of change of actin monomer concentration, Arp2/3 (Branch) conc., and filament barbed end
  // conc.
  CellScatra()->AddContributionToRHS(sources_);

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

  // remove sources from scatra rhs
  // this needs to be done because otherwise, we would always add neumann loads in every iteration
  sources_->Scale(-1.0);
  CellScatra()->AddContributionToRHS(sources_);
  sources_->Scale(-1.0);

  // set scalar transport values expressed in structural dofs
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // -------------------------------------------------------------------
  //                  evaluate new growth
  // -------------------------------------------------------------------
  // set displacement state
  StructureField()->Discretization()->SetState("displacement", StructureField()->Dispnp());

  // evaluate
  return EvaluateGrowth();
}


/*----------------------------------------------------------------------*
 |                                                          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::UpdateMatConf()
{
  if (myrank_ == 0) std::cout << "\n   Update Material Configuration ... " << std::endl;

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = AleToStructure(AleField()->Dispnp());

  // vector of current spatial displacements
  Teuchos::RCP<const Epetra_Vector> dispnp =
      StructureField()->Dispnp();  // change to ExtractDispn() for overlap

  // material displacements
  Teuchos::RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispnp->Map()), true);

  // set state
  StructureField()->Discretization()->SetState(0, "displacement", dispnp);

  // set state
  StructureField()->Discretization()->SetState(
      0, "material_displacement", StructAle()->GetMaterialDisplacementNpPtr());

  // calc difference between spatial structure deformation and ale deformation
  // disalenp = d_ale - d_struct
  disalenp->Update(-1.0, *dispnp, 1.0);

  double normdisale = -1234.0;
  disalenp->Norm2(&normdisale);

  // save the result in delta_ale_
  delta_ale_->Update(1.0, *disalenp, 0.0);

  // loop over all row nodes to fill graph
  for (int k = 0; k < StructureField()->Discretization()->NumMyRowNodes(); ++k)
  {
    int gid = StructureField()->Discretization()->NodeRowMap()->GID(k);
    DRT::Node* node = StructureField()->Discretization()->gNode(gid);
    DRT::Element** ElementPtr = node->Elements();
    int numelement = node->NumElement();

    const int numdof = StructureField()->Discretization()->NumDof(0, node);

    // create Xmat for 3D problems
    double XMat[numdof];
    double XMesh[numdof];

    for (int dof = 0; dof < numdof; ++dof)
    {
      int dofgid = StructureField()->Discretization()->Dof(0, node, dof);
      int doflid = (dispnp->Map()).LID(dofgid);
      XMesh[dof] = node->X()[dof] + (*dispnp)[doflid] + (*disalenp)[doflid];
    }

    // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
    AdvectionMap(XMat, XMesh, ElementPtr, numelement, true);

    // store in dispmat
    for (int dof = 0; dof < numdof; ++dof)
    {
      int dofgid = StructureField()->Discretization()->Dof(0, node, dof);
      int doflid = (dispnp->Map()).LID(dofgid);
      (*dismat)[doflid] = XMat[dof] - node->X()[dof];
    }
  }  // end row node loop

  // apply material displacements to structural field
  // if advection map is not succesful --> use old xmat
  StructAle()->UpdateMaterialDisplacements(dismat);

  return;
}


/*----------------------------------------------------------------------*
 |                                                          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::UpdateSpatConf()
{
  if (myrank_ == 0) std::cout << "\n   Update Spatial Configuration ... " << std::endl;

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = AleToStructure(AleField()->Dispnp());

  // get structure dispnp vector
  Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();

  // update per absolute vector
  dispnp->Update(1.0, *disalenp, 0.0);

  // synchronize nox state and global state
  // update velocities and accelerations
  StructureField()->SetState(dispnp);

  return;
}


/*----------------------------------------------------------------------*/
// advection map assembly analogous to wear framework       rauch 01/16 |
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::AdvectionMap(double* Xtarget,  // out
    double* Xsource,                                                      // in
    DRT::Element** ElementPtr,                                            // in
    int numelements,                                                      // in
    bool spatialtomaterial)                                               // in
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // define source and target configuration
  std::string sourceconf;
  std::string targetconf;

  if (spatialtomaterial)
  {
    sourceconf = "displacement";
    targetconf = "material_displacement";
  }
  else
  {
    sourceconf = "material_displacement";
    targetconf = "displacement";
  }

  // get state
  Teuchos::RCP<const Epetra_Vector> dispsource =
      StructureField()->Discretization()->GetState(sourceconf);
  Teuchos::RCP<const Epetra_Vector> disptarget =
      StructureField()->Discretization()->GetState(targetconf);

  // found element the spatial coordinate lies in
  bool found = false;

  // parameter space coordinates
  double e[3];
  double ge1 = 1e12;
  double ge2 = 1e12;
  double ge3 = 1e12;
  int gele = 0;

  // loop over adjacent elements
  for (int jele = 0; jele < numelements; jele++)
  {
    // get element
    DRT::Element* actele = ElementPtr[jele];

    // get element location vector, dirichlet flags and ownerships
    DRT::Element::LocationArray la(1);
    actele->LocationVector(*StructureField()->Discretization(), la, false);

    if (ndim == 2)
    {
      if (actele->Shape() == DRT::Element::quad4)
        WEAR::UTILS::av<DRT::Element::quad4>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::quad8)
        WEAR::UTILS::av<DRT::Element::quad8>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::quad9)
        WEAR::UTILS::av<DRT::Element::quad9>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::tri3)
        WEAR::UTILS::av<DRT::Element::tri3>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::tri6)
        WEAR::UTILS::av<DRT::Element::tri6>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else
        dserror("ERROR: shape function not supported!");

      // if parameter space coord. 'e' does not lie within any element (i.e. found = false),
      // then jele is the element lying closest near the considered spatial point.
      if (found == false)
      {
        if (abs(ge1) > 1.0 and abs(e[0]) < abs(ge1))
        {
          ge1 = e[0];
          gele = jele;
        }
        if (abs(ge2) > 1.0 and abs(e[1]) < abs(ge2))
        {
          ge2 = e[1];
          gele = jele;
        }
      }
    }
    else
    {
      if (actele->ElementType() == DRT::ELEMENTS::So_hex8Type::Instance() or
          actele->ElementType() == DRT::ELEMENTS::So_hex8ScatraType::Instance())
        WEAR::UTILS::av<DRT::Element::hex8>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
        WEAR::UTILS::av<DRT::Element::hex20>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
        WEAR::UTILS::av<DRT::Element::hex27>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
        WEAR::UTILS::av<DRT::Element::tet4>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
        WEAR::UTILS::av<DRT::Element::tet10>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else
        dserror("ERROR: element type not supported!");

      // if parameter space coord. 'e' does not lie within any element (i.e. found = false),
      // then jele is the element lying closest near the considered spatial point.
      if (found == false)
      {
        if (abs(ge1) > 1.0 and abs(e[0]) < abs(ge1))
        {
          ge1 = e[0];
          gele = jele;
        }

        if (abs(ge2) > 1.0 and abs(e[1]) < abs(ge2))
        {
          ge2 = e[1];
          gele = jele;
        }
        if (abs(ge3) > 1.0 and abs(e[2]) < abs(ge3))
        {
          ge3 = e[2];
          gele = jele;
        }
      }
    }

    // leave when element is found
    if (found == true) return;

  }  // end loop over adj elements

  // ****************************************
  //  if not displaced into elements: get
  //  Xtarget from closest element 'gele'
  // ****************************************
  DRT::Element* actele = ElementPtr[gele];

  // get element location vector, dirichlet flags and ownerships
  DRT::Element::LocationArray la(1);
  actele->LocationVector(*StructureField()->Discretization(), la, false);

  if (ndim == 2)
  {
    if (actele->Shape() == DRT::Element::quad4)
      WEAR::UTILS::av<DRT::Element::quad4>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::quad8)
      WEAR::UTILS::av<DRT::Element::quad8>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::quad9)
      WEAR::UTILS::av<DRT::Element::quad9>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::tri3)
      WEAR::UTILS::av<DRT::Element::tri3>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::tri6)
      WEAR::UTILS::av<DRT::Element::tri6>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else
      dserror("ERROR: shape function not supported!");
  }
  else
  {
    if (actele->ElementType() == DRT::ELEMENTS::So_hex8Type::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_hex8ScatraType::Instance())
      WEAR::UTILS::av<DRT::Element::hex8>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
      WEAR::UTILS::av<DRT::Element::hex20>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
      WEAR::UTILS::av<DRT::Element::hex27>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
      WEAR::UTILS::av<DRT::Element::tet4>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
      WEAR::UTILS::av<DRT::Element::tet10>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else
      dserror("ERROR: element type not supported!");
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | transform from structure to ale map                      rauch 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC_PROTRUSIONFORMATION::StructureToAle(
    Teuchos::RCP<const Epetra_Vector> vec)
{
  if (AleStruCoupling() == Teuchos::null) dserror("RCP to Coupling object points to Teuchos::null");
  return AleStruCoupling()->SlaveToMaster(vec);
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      rauch 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC_PROTRUSIONFORMATION::AleToStructure(
    Teuchos::RCP<const Epetra_Vector> vec)
{
  return AleStruCoupling()->MasterToSlave(vec);
}


/*----------------------------------------------------------------------*
 | Calculate source/sink values                             rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::EvaluateSources()
{
  // get parameter list
  Teuchos::ParameterList params;

  // get time step
  int curr_step = scatra_->ScaTraField()->Step();

  // add time step
  params.set<int>("current step", curr_step);
  // add scatra discretization
  params.set<std::string>("scatradisname", "cellscatra");
  // add number of disp dofset
  params.set<int>("ndsdisp", 1);
  // add time step size
  params.set<double>("dt", Dt());

  // set states
  scatra_->ScaTraField()->Strategy()->SetState(0, "phinp", scatra_->ScaTraField()->Phinp());
  scatra_->ScaTraField()->Strategy()->SetState(0, "phin", scatra_->ScaTraField()->Phin());

  // add action for growth evaluation
  params.set<int>("action", SCATRA::calc_cell_growth_sourcesandsinks);
  // evaluate condition on auxiliary discretization in strategy
  scatra_->ScaTraField()->Strategy()->EvaluateCondition(params, Teuchos::null, Teuchos::null,
      sources_, Teuchos::null, Teuchos::null, "ScatraHeteroReactionSlave", -1);

  // initialize norm
  double sourcenorm = -1234.0;
  // calc norm
  sources_->Norm2(&sourcenorm);
  // print information
  std::cout << "\n   Norm of growth related source vector " << std::setprecision(10) << sourcenorm
            << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 | Calculate growth values from source/sink values          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::EvaluateGrowth()
{
  //////////////////////////////////////////////////////////////////////////////////////
  // Evaluate biochemical transport and reaction at growth surface.
  // Evaluate Condition (FSICoupling).
  // Evaluate polymerisation (change of monomer conc. and change of pointed end conc.
  // at conditioned surface) and evaluate growth due to polymerisation.
  //////////////////////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  // add action for growth evaluation
  params.set<std::string>("action", "calc_cell_growth");
  // add timestep to parameterlit
  params.set<double>("dt", Dt());

  if (not conditiondofrowmap_->UniqueGIDs())
    dserror("conditiondofrowmap_ is not unique! Something went wrong!");

  // least squares system-matrix
  Teuchos::RCP<LINALG::SparseOperator> leastsquares_matrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*conditiondofrowmap_, 18, false, true));
  // least squares right-hand side
  Teuchos::RCP<Epetra_Vector> leastsquares_rhs = LINALG::CreateVector(*conditiondofrowmap_, true);
  // least squares error optimal nodal growth vector
  Teuchos::RCP<Epetra_Vector> leastsquares_growth =
      LINALG::CreateVector(*conditiondofrowmap_, true);

  // set states
  StructureField()->Discretization()->SetState(1, "phinp", scatra_->ScaTraField()->Phinp());
  StructureField()->Discretization()->SetState(1, "phin", scatra_->ScaTraField()->Phin());
  StructureField()->Discretization()->SetState(1, "rates", sources_);

  // evaluate least squares growth todo remove auxiliary use of FSICoupling condition
  StructureField()->Discretization()->EvaluateCondition(params, leastsquares_matrix, Teuchos::null,
      leastsquares_rhs, Teuchos::null, Teuchos::null, "FSICoupling", -1);

  /////////////////////////////////////////////////////
  // Calc nodal growth from gauss point growth with
  // least squares method, i.e.:
  // Minimize sum|N * d - d_gp|²
  // And solve the resulting system of equations for
  // the leastsquares_growth x :
  // b = A * x --> x = A^1 * b
  // with
  // b = leastsquares_rhs from Evaluate Condition
  // A = leastsquares_matrix from Evaluate Condition
  ///////////////////////////////////////////////////
  // SOLVE LEAST SQUARES SYSTEM
  // solver setup
  Teuchos::ParameterList param_solve = DRT::Problem::Instance()->UMFPACKSolverParams();
  bool refactor = true;
  bool reset = true;
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(
      new LINALG::Solver(param_solve, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(solver->Params());

  // complete matrix
  leastsquares_matrix->Complete();

  // solve for least squares optimal nodal values
  solver->Solve(leastsquares_matrix->EpetraOperator(), leastsquares_growth, leastsquares_rhs,
      refactor, reset);

  // find growth values on growth surface
  int numdof = conditiondofrowmap_->NumMyElements();
  int numdofgrowth = growth_step_->MyLength();
  if (numdof != numdofgrowth) dserror("size of growth surface not matching");

  for (int i = 0; i < numdofgrowth; i++)
  {
    // fill step increment (growth from t_n -> t_{n+1})
    growth_step_->ReplaceMyValue(i, 0, leastsquares_growth->Values()[i]);
  }

  // update new growth
  int err = growth_np_->Update(1.0, *growth_n_, 1.0, *growth_step_, 0.0);
  if (err != 0) dserror("epetra update of  growth_np_ returned err=%d", err);

  return;
}


/*----------------------------------------------------------------------*
 | convergence check                                        rauch 01/16 |
 *----------------------------------------------------------------------*/
bool SSI::SSI_Part2WC_PROTRUSIONFORMATION::ConvergenceCheck(int itnum)
{
  // initialize flags
  bool stopnonliniter = false;
  bool SSI_converged = false;

  // norm of growth
  double growthnorm_n = -1234.0;
  double growthincnorm = -1234.0;
  double growthnorm = -1234.0;

  // do the base class convergence check
  SSI_converged = SSI::SSI_Part2WC::ConvergenceCheck(itnum);

  // set converged if base class check returned true
  if (SSI_converged) stopnonliniter = true;

  // get the new the growth increment
  //  growthinc_->Update(-1.0,*growth_n_,1.0);
  //  growthinc_->Update(-1.0,*growth_step_,1.0);
  int err = growthinc_->Update(1.0, *growth_np_, -1.0, *x_, 0.0);
  if (err != 0) dserror("epetra update returned err=%d", err);

  // update state
  err = x_->Update(omega_, *growthinc_, 1.0);
  if (err != 0) dserror("epetra update returned err=%d", err);

  // get L2-Norm of growth
  growthinc_->Norm2(&growthincnorm);

  // get L2-Norm of growth
  growth_step_->Norm2(&growthnorm);

  // get L2-Norm of growth at old time step
  x_->Norm2(&growthnorm_n);

  // proc 0 prints relative growth increment
  if (myrank_ == 0)
  {
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
    printf("L2-Norm of state                  =   %10.6E \n", growthnorm_n);
    printf("L2-Norm of growth from tn->tn+1   =   %10.6E \n", growthnorm);
    printf("L2-Norm of growt increment i->i+1 =  %10.3E < %10.3E",
        growthincnorm / growthinc_->GlobalLength(), ittol_);
    printf("\n");
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
  }

  // ignore base class check if growth is converged
  if (DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE"),
          "SSICOUPVARIABLE") == INPAR::CELL::coup_growth_growth)
  {
    if ((growthincnorm / growthinc_->GlobalLength()) <= ittol_)
    {
      if (myrank_ == 0) std::cout << "Growth is converged. Ignore convergence of SSI." << std::endl;

      stopnonliniter = true;
    }
  }
  else if (DRT::INPUT::IntegralValue<int>(
               DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE"),
               "SSICOUPVARIABLE") == INPAR::CELL::coup_growth_undefined)
    dserror(
        "set SSICOUPVARIABLE in section ---CELL DYNAMIC/PROTRUSION MODULE to 'growth' or 'ssi'");

  // tell if we converged
  return stopnonliniter;
}


/*------------------------------------------------------------------------*
 | BuildConditionDofRowMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::BuildConditionDofRowMap(const DRT::Condition* condition,
    const Teuchos::RCP<const DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& conddofmap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();
  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < (int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));

    // if node is owned by calling proc
    if (dis->NodeRowMap()->LID(currnode->Id()) != -1)
    {
      std::vector<int> dofs = dis->Dof(0, currnode);

      for (int dim = 0; dim < 3; ++dim)
      {
        mydirichdofs.push_back(dofs[dim]);
      }
    }  // if node owned by calling proc
  }    // loop over all conditioned nodes

  int nummydirichvals = mydirichdofs.size();
  conddofmap =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}


/*------------------------------------------------------------------------*
 | update the current states in every iteration               rauch 05/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::IterUpdateStates()
{
  // perform the update from the base class
  SSI::SSI_Part2WC::IterUpdateStates();

  // clear sources vector
  sources_->PutScalar(0.0);

  // clear growth step vector
  growth_step_->PutScalar(0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::UpdateAndOutput()
{
  // call update and output from base class
  SSI::SSI_Part2WC::UpdateAndOutput();

  // update growth_{n}
  int err = growth_n_->Update(1.0, *x_, 0.0);
  if (err != 0) dserror("update of growth_n_ returned err=%d", err);

  // clear growth step increment
  growth_step_->PutScalar(0.0);

  // clear sources vector
  sources_->PutScalar(0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::SetScatraWrapper(
    Teuchos::RCP<::ADAPTER::AdapterScatraWrapperCellMigration> scatra_wrapper)
{
  cell_scatra_wrapper_ = scatra_wrapper;
}
