/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.cpp

\brief specialization of ssi2wc, including "structale"-surface growth

\level 3

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

 *----------------------------------------------------------------------*/
#include "ssi_partitioned_2wc_protrusionformation.H"
#include "../drt_immersed_problem/immersed_field_exchange_manager.H"
#include "../drt_inpar/inpar_cell.H"

#include "../drt_adapter/ad_ale_fluid.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_nodematchingoctree.H"

#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex20.H"
#include "../drt_so3/so_hex27.H"
#include "../drt_so3/so_tet4.H"
#include "../drt_so3/so_tet10.H"
#include "../drt_so3/so3_scatra.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_io/io_control.H"
#include "../drt_wear/wear_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_clonestrategy.H"


/*----------------------------------------------------------------------*
 | constructor                                              rauch 01/16 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_PROTRUSIONFORMATION::SSI_Part2WC_PROTRUSIONFORMATION(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
: SSI_Part2WC(comm, globaltimeparams, scatraparams, structparams,struct_disname,scatra_disname)
{
  // set communicator
  myrank_ = comm.MyPID();

  // check if scatra in cell is set up with ale description
  if(not ScaTraField()->ScaTraField()->IsALE())
    dserror("We need an ALE description for the cell-scatra field!");

  // additional setup for ale
  SetupDiscretizations(comm,"cell","cellscatra");

  // initialize pointer to structure
  specialized_structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(StructureField());
  if(specialized_structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(new ADAPTER::AleBaseAlgorithm(DRT::Problem::Instance()->SSIControlParams(), DRT::Problem::Instance()->GetDis("ale")));
  ale_ =  Teuchos::rcp_dynamic_cast<ADAPTER::AleFluidWrapper>(ale->AleField());
  if(ale_ == Teuchos::null)
    dserror("cast from ADAPTER::Ale to ADAPTER::AleFsiWrapper failed");
  // create empty operator
  AleField()->CreateSystemMatrix();

  // build coupling objects for dof transfer between structure and ale
  const int ndim = DRT::Problem::Instance()->NDim();

  // create ale-struct coupling
  const Epetra_Map* celldofmap =
      StructureField()->Discretization()->NodeRowMap();
  const Epetra_Map* aledofmap = AleField()->Discretization()->NodeRowMap();

  // if there are two identical nodes (i.e. for initial contact) the nodes matching creates an error !!!
  coupalestru_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupalestru_->SetupCoupling(*AleField()->Discretization(),
      *StructureField()->Discretization(), *aledofmap, *celldofmap, ndim, true, 1e-06);

  // set-up aux FSI interface
  Teuchos::RCP<STR::AUX::MapExtractor> interface = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface->Setup(*(scatra_->ScaTraField()->Discretization()), *(scatra_->ScaTraField()->Discretization()->DofRowMap()) );

  // growth state
  x_=Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth from last time step
  growth_n_=Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth at new time step
  growth_np_=Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // growth step increment growth_{n+1} = growth_{n} + growth_{step}
  growth_step_=Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true);

  // initialize growth increment vector
  growthinc_=Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(AleField()->Interface()->cond_fsi)), true),

  // initialize delta_ale_ vector
  delta_ale_=Teuchos::rcp(new Epetra_Vector(AleField()->Dispnp()->Map(), true));

  // source/sink vector to be added to scatra rhs
  sources_ = LINALG::CreateVector(*(scatra_->ScaTraField()->Discretization()->DofRowMap()),true);

  // numdof actin
  numdof_actin_ = DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE").get<int>("NUMDOF_ACTIN");
  if(myrank_==0)
    std::cout<<"\n  Number of Actin Dof in Scalar Transport System: "<<numdof_actin_<<"\n"<<std::endl;

  // get pointer to the ImmersedFieldExchangeManager
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

  // set pointer to the concentrations at n+1
  exchange_manager_->SetPointerToPhinps(scatra_->ScaTraField()->Phinp());

  // initialize pointer to phin at n
  Teuchos::RCP<Epetra_MultiVector> phin = scatra_->ScaTraField()->Phin();
  exchange_manager_->SetPointerToPhins(phin);

  // set pointers to multivectors for the rates
  exchange_manager_->SetPointerToRates(sources_); // rates of polymerized actin monomers, barbed ends, and branches

  // mapping from struct bdry ele gids to scatra bdry ele gids
  structscatraelemap_=Teuchos::rcp(new std::map<int,int> );

  // relaxation
  omega_ = DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE").get<double>("RELAX_GROWTH");
  if(myrank_==0)
    std::cout<<"\n  Use fixed relaxation parameter for protrusion formation omega="<<omega_<<"\n"<<std::endl;

  // NODE MATCHING
  // do the matching
  MatchNodes(*(StructureField()->Discretization()),
             *(scatra_->ScaTraField()->Discretization()),
             structnodemap_,
             scatranodemap_,
             "CellSurfVolCoupling");

  // build an element map associating matched elements of struct and scatra discretisations
  BuildMasterGeometryToSlaveDisEleMap(*(StructureField()->Discretization()),
                                      *(scatra_->ScaTraField()->Discretization()),
                                      structnodemap_,
                                      scatranodemap_,
                                      structscatraelemap_,
                                      "CellSurfVolCoupling");
  // build map
  BuildConditionDofRowMap((StructureField()->Discretization())->GetCondition("FSICoupling"),StructureField()->Discretization(),conditiondofrowmap_);
  // print map
  std::cout<<"MAPPING CELL 'CellSurfVolCoupling' GEOMETRY TO SCATRA 'CellSurfVolCoupling' ELEMENTS ... "<<std::endl;
  std::cout<<structscatraelemap_->size()<<" MAPPED ELEMENTS.\n"<<std::endl;
  std::cout<<"STRUCT ELE IDs --> SCATRA ELE IDs"<<std::endl;
  for(std::map<int,int>::iterator it=structscatraelemap_->begin(); it!=structscatraelemap_->end();++it)
    std::cout<<"    "<<it->first<<"               "<<it->second<<std::endl;
}


/*----------------------------------------------------------------------*
 | Solve structure filed                                    rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
    << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
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
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());

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
  double normofgrowth=-1234.0;
  // calc norm
  growthincrement->Norm2(&normofgrowth);

  if(myrank_==0)
    std::cout<<"=======================================\n"
               "DoAleStep\n"
               "Norm of applied growth "<<std::setprecision(10)<<normofgrowth<<"\n"<<
               "======================================="<<std::endl;

  // get lagrangian structure displacements
  Teuchos::RCP<Epetra_Vector> dispnpstru = StructureToAle(StructureField()->Dispnp());

  // update ale field with lagrangian structure displacements
  AleField()->WriteAccessDispnp()->Update(1.0, *dispnpstru, 0.0);

  // application of interface displacements as dirichlet conditions
  AleField()->ApplyInterfaceDisplacements(growthincrement);

  // solve time step
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
    std::cout
        << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // rate of change of actin monomer concentration, Arp2/3 (Branch) conc., and filament barbed end conc.
  scatra_->ScaTraField()->AddContributionToRHS(sources_);

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

  // remove sources from scatra rhs
  // this needs to be done because otherwise, we would always add neumann loads in every iteration
  sources_->Scale(-1.0);
  scatra_->ScaTraField()->AddContributionToRHS(sources_);
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
  if(myrank_==0)
    std::cout<<"\n   Update Material Configuration ... "<<std::endl;

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = AleToStructure(AleField()->Dispnp());

  // vector of current spatial displacements
  Teuchos::RCP<const Epetra_Vector> dispnp = StructureField()->Dispnp(); // change to ExtractDispn() for overlap

  // material displacements
  Teuchos::RCP<Epetra_Vector> dismat = Teuchos::rcp(
      new Epetra_Vector(dispnp->Map()), true);

  // set state
  StructureField()->Discretization()->SetState(0, "displacement", dispnp);

  // set state
  StructureField()->Discretization()->SetState(0, "material_displacement",StructureField()->DispMat());

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

    const int numdof = StructureField()->Discretization()->NumDof(node);

    // create Xmat for 3D problems
    double XMat[numdof];
    double XMesh[numdof];

    for(int dof = 0; dof < numdof; ++dof)
    {
      int dofgid = StructureField()->Discretization()->Dof(node,dof);
      int doflid = (dispnp->Map()).LID(dofgid);
      XMesh[dof] = node->X()[dof] + (*dispnp)[doflid] + (*disalenp)[doflid];
    }

    // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
    AdvectionMap(XMat, XMesh, ElementPtr, numelement, true);

    // store in dispmat
    for(int dof = 0; dof < numdof; ++dof)
    {
      int dofgid = StructureField()->Discretization()->Dof(node,dof);
      int doflid = (dispnp->Map()).LID(dofgid);
      (*dismat)[doflid] = XMat[dof] - node->X()[dof];
    }
  } // end row node loop

  // apply material displacements to structural field
  // if advection map is not succesful --> use old xmat
  StructureField()->ApplyDisMat(dismat);

  return;
}


/*----------------------------------------------------------------------*
 |                                                          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::UpdateSpatConf()
{
  if(myrank_==0)
    std::cout<<"\n   Update Spatial Configuration ... "<<std::endl;

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = AleToStructure(AleField()->Dispnp());

  // get structure dispnp vector
  Teuchos::RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();

  // update per absolute vector
  dispnp->Update(1.0, *disalenp, 0.0);

  // also update velocity to be consistent with new dispnp
  // time step size
  const double dt = Dt();
  // OST Parameter
  const double theta = DRT::Problem::Instance()->CellMigrationParams().sublist("STRUCTURAL DYNAMIC").sublist("ONESTEPTHETA").get<double>("THETA");
  if(theta!=1.0)
    dserror("algorithm has not been tested with theta != 1.0");

  // new end-point velocities
  StructureField()->WriteAccessVelnp()->Update(1.0/(theta*dt), *dispnp,
                                              -1.0/(theta*dt), *((StructureField()->Dispn())),
                                               0.0);
  StructureField()->WriteAccessVelnp()->Update(-(1.0-theta)/theta, *((StructureField()->Veln())),
                                               1.0);
  return;
}


/*----------------------------------------------------------------------*/
// advection map assembly analogous to wear framework       rauch 01/16 |
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::AdvectionMap(
    double* Xtarget,            // out
    double* Xsource,            // in
    DRT::Element** ElementPtr,  // in
    int numelements,            // in
    bool spatialtomaterial)     // in
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // define source and target configuration
  std::string sourceconf;
  std::string targetconf;

  if(spatialtomaterial)
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
  Teuchos::RCP<const Epetra_Vector> dispsource = StructureField()->Discretization()->GetState(sourceconf);
  Teuchos::RCP<const Epetra_Vector> disptarget = StructureField()->Discretization()->GetState(targetconf);

  // found element the spatial coordinate lies in
  bool found = false;

  // parameter space coordinates
  double e[3];
  double ge1 = 1e12;
  double ge2 = 1e12;
  double ge3 = 1e12;
  int gele   = 0;

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
        WEAR::UTILS::av<DRT::Element::quad4>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::quad8)
        WEAR::UTILS::av<DRT::Element::quad8>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::quad9)
        WEAR::UTILS::av<DRT::Element::quad9>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::tri3)
        WEAR::UTILS::av<DRT::Element::tri3>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->Shape() == DRT::Element::tri6)
        WEAR::UTILS::av<DRT::Element::tri6>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
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
        WEAR::UTILS::av<DRT::Element::hex8>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
        WEAR::UTILS::av<DRT::Element::hex20>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
        WEAR::UTILS::av<DRT::Element::hex27>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
        WEAR::UTILS::av<DRT::Element::tet4>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
        WEAR::UTILS::av<DRT::Element::tet10>(actele, Xtarget, Xsource, dispsource, disptarget,
            la[0].lm_, found, e);
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
    if (found == true)
      return;

  } // end loop over adj elements

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
      WEAR::UTILS::av<DRT::Element::quad4>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::quad8)
      WEAR::UTILS::av<DRT::Element::quad8>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::quad9)
      WEAR::UTILS::av<DRT::Element::quad9>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::tri3)
      WEAR::UTILS::av<DRT::Element::tri3>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->Shape() == DRT::Element::tri6)
      WEAR::UTILS::av<DRT::Element::tri6>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else
      dserror("ERROR: shape function not supported!");
  }
  else
  {
    if (actele->ElementType() == DRT::ELEMENTS::So_hex8Type::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_hex8ScatraType::Instance())
      WEAR::UTILS::av<DRT::Element::hex8>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
      WEAR::UTILS::av<DRT::Element::hex20>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
      WEAR::UTILS::av<DRT::Element::hex27>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
      WEAR::UTILS::av<DRT::Element::tet4>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
      WEAR::UTILS::av<DRT::Element::tet10>(actele, Xtarget, Xsource, dispsource, disptarget,
          la[0].lm_, found, e);
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
    Teuchos::RCP<Epetra_Vector> vec)
{
  return AleStruCoupling()->SlaveToMaster(vec);
}


/*----------------------------------------------------------------------*
 | transform from structure to ale map                      rauch 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC_PROTRUSIONFORMATION::StructureToAle(
    Teuchos::RCP<const Epetra_Vector> vec)
{
  if (AleStruCoupling()==Teuchos::null)
    dserror("RCP to Coupling object points to Teuchos::null");
  return AleStruCoupling()->SlaveToMaster(vec);
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      rauch 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC_PROTRUSIONFORMATION::AleToStructure(
    Teuchos::RCP<Epetra_Vector> vec)
{
  return AleStruCoupling()->MasterToSlave(vec);
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
 |                                                          rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string scatra_disname)
{
  // call SetupDiscretizations in base class
  // Done in constructor of ssi_base

  // new ale part
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);

  // clone ale from structure for ssi-actin assembly
  // access the ale discretization
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = DRT::Problem::Instance()->GetDis("ale");
  if (!aledis->Filled()) aledis->FillComplete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->NumGlobalNodes()==0)
    dserror("ERROR: Structure discretization is empty!");

  // clone ale mesh from structure discretization
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(structdis,aledis);
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else
    dserror("ERROR: Reading an ALE mesh from the input file is not supported for this problem type.");

}


/*----------------------------------------------------------------------*
 | Calculate source values                                  rauch 01/16 |
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
  params.set<int>("ndsdisp",1);
  // add time step size
  params.set<double>("dt",Dt());

  // add action for growth evaluation
  params.set<int>("action",SCATRA::calc_cell_growth);
  // evaluate condition
  scatra_->ScaTraField()->Discretization()->EvaluateCondition(params,
                                                              Teuchos::null,
                                                              Teuchos::null,
                                                              sources_,
                                                              Teuchos::null,
                                                              Teuchos::null,
                                                              "CellSurfVolCoupling",
                                                              -1);

  // initialize norm
  double sourcenorm=-1234.0;
  // calc norm
  sources_->Norm2(&sourcenorm);
  // print information
  std::cout<<"\n   Norm of growth related source vector "<<std::setprecision(10)<<sourcenorm<<std::endl;

  return;
}


/*----------------------------------------------------------------------*
 | Calculate growth values                                  rauch 01/16 |
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

  // add master (structure) nodemap to parameterlist
  params.set<Teuchos::RCP<Epetra_Map> >("masternodemap",structnodemap_);
  // add slave (scatra) nodemap to parameterlist
  params.set<Teuchos::RCP<Epetra_Map> >("slavenodemap",scatranodemap_);
  // add action for growth evaluation
  params.set<std::string>("action","calc_cell_growth");
  // add scatra discretization
  params.set<std::string>("scatradisname", "cellscatra");
  // add condition string
  params.set<std::string>("condstring", "FSICoupling");
  // add rcp to struct<->scatra brdy ele map to parameterlist
  params.set<Teuchos::RCP<std::map<int,int> > >("structscatraelemap",structscatraelemap_);
  // add timestep to parameterlit
  params.set<double>("dt",Dt());

  // least squares system-matrix
  Teuchos::RCP<LINALG::SparseOperator> leastsquares_matrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*conditiondofrowmap_,18,false,true));
  // least squares right-hand side
  Teuchos::RCP<Epetra_Vector> leastsquares_rhs = LINALG::CreateVector(*conditiondofrowmap_,true);
  // least squares error optimal nodal growth vector
  Teuchos::RCP<Epetra_Vector> leastsquares_growth = LINALG::CreateVector(*conditiondofrowmap_,true);

  // evaluate least squares growth
  StructureField()->Discretization()->EvaluateCondition(params,
      leastsquares_matrix,
      Teuchos::null,
      leastsquares_rhs,
      Teuchos::null,
      Teuchos::null,
      "CellSurfVolCoupling",
      -1);

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
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(param_solve, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(solver->Params());

  // complete matrix
  leastsquares_matrix->Complete();

  // solve for least squares optimal nodal values
  solver->Solve(leastsquares_matrix->EpetraOperator(), leastsquares_growth, leastsquares_rhs, refactor, reset);

  // find growth values on growth surface
  int numdof = conditiondofrowmap_->NumGlobalElements();
  int numdofgrowth = growth_step_->MyLength();
  if (numdof != numdofgrowth)
    dserror("size of growth surface not matching");

  for (int i=0; i<numdofgrowth; i++)
  {
    // fill step increment (growth from t_n -> t_{n+1})
    growth_step_->ReplaceMyValue(i, 0, leastsquares_growth->Values()[i]);
  }

  // update new growth
  int err = growth_np_->Update(1.0,*growth_n_,1.0,*growth_step_,0.0);
  if(err!=0)
    dserror("epetra update of  growth_np_ returned err=%d",err);

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
  if(SSI_converged)
    stopnonliniter=true;

  // get the new the growth increment
//  growthinc_->Update(-1.0,*growth_n_,1.0);
//  growthinc_->Update(-1.0,*growth_step_,1.0);
  int err = growthinc_->Update(1.0,*growth_np_,-1.0,*x_,0.0);
  if(err!=0)
    dserror("epetra update returned err=%d",err);

  // update state
  err = x_->Update(omega_,*growthinc_,1.0);
  if(err!=0)
    dserror("epetra update returned err=%d",err);

  // get L2-Norm of growth
  growthinc_->Norm2(&growthincnorm);

  // get L2-Norm of growth
  growth_step_->Norm2(&growthnorm);

  // get L2-Norm of growth at old time step
  x_->Norm2(&growthnorm_n);

  // proc 0 prints relative growth increment
  if(myrank_==0)
  {
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
    printf("L2-Norm of state                  =   %10.6E \n", growthnorm_n);
    printf("L2-Norm of growth from tn->tn+1   =   %10.6E \n", growthnorm);
    printf("L2-Norm of growt increment i->i+1 =  %10.3E < %10.3E", growthincnorm/growthinc_->GlobalLength(), ittol_);
    printf("\n");
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
  }

  // ignore base class check if growth is converged
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE"),"COUPVARIABLE") == INPAR::CELL::coup_growth_growth)
  {
    if((growthincnorm/growthinc_->GlobalLength())<=ittol_)
    {
      if(myrank_==0)
        std::cout<<"Growth is converged. Ignore convergence of SSI."<<std::endl;

      stopnonliniter=true;
    }
  }
  else if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CellMigrationParams().sublist("PROTRUSION MODULE"),"COUPVARIABLE") == INPAR::CELL::coup_growth_undefined)
    dserror("set COUPVARIABLE in section ---CELL DYNAMIC/PROTRUSION MODULE to 'growth' or 'ssi'");

  // tell if we converged
  return stopnonliniter;
}


/*------------------------------------------------------------------------*
 | matching nodes of master and slave field                   rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::MatchNodes(
    DRT::Discretization& masterdis,
    DRT::Discretization& slavedis,
    Teuchos::RCP<Epetra_Map>& masternodemap_matched,
    Teuchos::RCP<Epetra_Map>& slavenodemap_matched,
    const std::string& condname)
{
  // get conditioned master nodes
  const std::vector<int>* masternodes = masterdis.GetCondition(condname)->Nodes();

  // get conditioned slave nodes
  const std::vector<int>* slavenodes = slavedis.GetCondition(condname)->Nodes();

  // sanity check
  if(masternodes->size()!=slavenodes->size())
    dserror("current algorithm relies on matching nodes of scatra surface discretization and cell boundary discretization!");

  // node vector of matched nodes
  std::vector<int> permslavenodes;

  // match master and slave nodes using Peter's octtree
  DRT::UTILS::NodeMatchingOctree tree(masterdis, *masternodes, 150, 1e-8);

  // coupled nodes - find matching nodes between master and slave dis.
  // coupled = nodes are at same location
  std::map<int,std::pair<int,double> > coupling;
  tree.FindMatch(slavedis, *slavenodes, coupling);

  // extract permutation
  std::vector<int> patchedmasternodes;
  patchedmasternodes.reserve(coupling.size());
  permslavenodes.reserve(slavenodes->size());

  for (unsigned i=0; i<masternodes->size(); ++i)
  {
    // get master node gid
    int gid = (*masternodes)[i];

    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour, the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int,double>& coupled = coupling[gid];
#if 0
      if (coupled.second > 1e-7)
        dserror("Coupled nodes (%d,%d) do not match. difference=%e", gid, coupled.first, coupled.second);
#endif
      patchedmasternodes.push_back(gid);       //< coupled masternodes
      permslavenodes.push_back(coupled.first); //< coupled slavenodes
    }
  }

  // epetra maps in original distribution
  // new masternode map (with coupled node IDs)
  masternodemap_matched = Teuchos::rcp(new Epetra_Map(
                     -1,
                     masternodes->size(),
                     &patchedmasternodes[0],
                     0,
                     masterdis.Comm()));

  // new slavenode map (with coupled node IDs)
  slavenodemap_matched = Teuchos::rcp(new Epetra_Map(
                    -1, permslavenodes.size(),
                    &permslavenodes[0],
                    0,
                    slavedis.Comm()));
  // new slavenodes
  std::vector<int> permslavenodes_new (slavenodemap_matched->MyGlobalElements(),
      slavenodemap_matched->MyGlobalElements() + slavenodemap_matched->NumMyElements());

} // MatchNodes


/*------------------------------------------------------------------------*
 | matching nodes of master and slave field                   rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::BuildMasterGeometryToSlaveDisEleMap(
    const DRT::Discretization& masterdis,
    const DRT::Discretization& slavedis,
    const Teuchos::RCP<Epetra_Map>& masternodemap,
    const Teuchos::RCP<Epetra_Map>& slavenodemap,
    Teuchos::RCP<std::map<int,int> > maptofill,
    const std::string& condname)
{
  if(maptofill==Teuchos::null)
    dserror("RCP to map object needs to be provided to this method");

  if(masternodemap->NumGlobalElements()==0)
    dserror("masternodemap is empty");

  if(slavenodemap->NumGlobalElements()==0)
    dserror("slavenodemap is empty");

  int structeleID=-1234; //< structure element ID to be matched -> key value
  int scatraeleID=-1234; //< scatra element ID to be matched -> second value

  // get geometry of condition 'condname' on 'masterdis'
  std::map<int,Teuchos::RCP<DRT::Element> >& mastergeom = masterdis.GetCondition(condname)->Geometry();

  // loop over master geometry elements
  for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator it=mastergeom.begin(); it!=mastergeom.end(); ++it)
  {
    // get current structure element ID
    structeleID = it->second->Id();
    // get number of nodes of current structure element
    int numnode = it->second->NumNode();
    if(numnode!=4)
      dserror("in this context, we expect only geometry elements with 4 nodes");

    // create vector
    std::vector<int> adjacentslaveeleIds_gathered;
    // prepare vector
    if(adjacentslaveeleIds_gathered.size())
      adjacentslaveeleIds_gathered.clear();

    // loop over master geometry element nodes
    for (int node=0; node<numnode; node++)
    {
      // find corresponding scatra node ID
      // gid of current structure element node
      int gid_master_node = (it->second->NodeIds())[node];
      // lid of current master node in masternodemap
      int lid_master_node = masternodemap->LID(gid_master_node);
      // gid of corresponding scatra element node (slavenodemap == nodemap of conditioned nodes - corresponds to currscatra elements)
      int gid_slave_node = slavenodemap->GID(lid_master_node);
      // get scatra node corresponding to structure node (with gidscatra)
      DRT::Node* node_slave = slavedis.gNode(gid_slave_node);

      ///////////////////////////////////////////////////
      // FIND SLAVE ELEMENT MATCHING MASTER ELEMENT
      ///////////////////////////////////////////////////
      // get pointer to adjacent elements of slave node
      DRT::Element** ElementPtr = node_slave->Elements();
      if (ElementPtr == NULL)
        dserror("could not get element pointer");

      // loop over all adjacent elements
      for (int jele=0; jele<node_slave->NumElement(); jele++)
      {
        // write IDs of adjacent elements in vector
        int adjacentslaveeleId = ElementPtr[jele]->Id();
        // we can only match elements with equal number of nodes
        if(ElementPtr[jele]->NumNode()==numnode)
        {
          adjacentslaveeleIds_gathered.push_back(adjacentslaveeleId);
        }
      }// end loop over all adjacent elements

    }//end loop over master geometry element nodes

    // FIND MATCHING ELEMENT ID OF ELEMENT ID VECTOR
    // count variable
    int count = 0;
    // if a slave element is an adjacent element to every node of the master element,
    // then it is the matching slave element (count should then be numnode-1)
    const int check = numnode-1;

    // loop over adjacent element ID vector
    for (unsigned int i=0; i<adjacentslaveeleIds_gathered.size(); i++)
    {

      for (unsigned int j=i+1; j<adjacentslaveeleIds_gathered.size(); j++)
      {
        if (adjacentslaveeleIds_gathered[i] == adjacentslaveeleIds_gathered[j])
        {
          count ++;
        }
      }

      if (count == check)
      {
        scatraeleID = adjacentslaveeleIds_gathered[i];
        break;
      }
      count = 0;
    }// end loop over element ID vector

    // fill element map (first : ID of structure element, second: ID of corresponding scatra element
    maptofill->insert( std::pair<int,int>(structeleID, scatraeleID) );

  }// end loop over master geometry elements

  return;
}// end BuildStructToScatraEleMap


/*------------------------------------------------------------------------*
 | BuildConditionDofRowMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_PROTRUSIONFORMATION::BuildConditionDofRowMap(const DRT::Condition* condition,
                                                                   const Teuchos::RCP<const DRT::Discretization> dis,
                                                                   Teuchos::RCP<Epetra_Map>& conddofmap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<(int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));
    std::vector<int> dofs = dis->Dof(0,currnode);

    for (int dim=0;dim<3;++dim)
    {
      mydirichdofs.push_back(dofs[dim]);
    }
  }

  int nummydirichvals = mydirichdofs.size();
  conddofmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

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
  int err = growth_n_->Update(1.0,*x_,0.0);
  if(err!=0)
    dserror("update of growth_n_ returned err=%d",err);

  // clear growth step increment
  growth_step_->PutScalar(0.0);

  // clear sources vector
  sources_->PutScalar(0.0);

  return;
}
