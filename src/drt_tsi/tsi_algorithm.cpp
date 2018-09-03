/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.cpp

\brief  Basis of all TSI algorithms that perform a coupling between the linear
        momentum equation and the heat conduction equation
\level 2
<pre>
  \maintainer  Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_algorithm.H"
#include "tsi_defines.H"
#include "tsi_utils.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"

#include "../drt_adapter/adapter_thermo.H"
#include "../drt_lib/drt_discret.H"

// for coupling of nonmatching meshes
#include "../drt_adapter/adapter_coupling_volmortar.H"
#include "../drt_volmortar/volmortar_utils.H"

// contact
#include "../drt_contact/contact_lagrange_strategy.H"
#include "../drt_contact/contact_tsi_lagrange_strategy.H"
#include "../drt_contact/contact_nitsche_strategy_tsi.H"
#include "../drt_contact/meshtying_contact_bridge.H"
#include "../drt_contact/contact_strategy_factory.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"
#include "../drt_structure_new/str_model_evaluator_structure.H"
#include "../drt_mortar/mortar_multifield_coupling.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(const Epetra_Comm& comm)
    : AlgorithmBase(comm, DRT::Problem::Instance()->TSIDynamicParams()),
      dispnp_(Teuchos::null),
      tempnp_(Teuchos::null),
      matchinggrid_(DRT::INPUT::IntegralValue<bool>(
          DRT::Problem::Instance()->TSIDynamicParams(), "MATCHINGGRID")),
      volcoupl_(Teuchos::null)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermodis = DRT::Problem::Instance()->GetDis("thermo");

  // get the problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  // get the restart step
  const int restart = problem->Restart();

  if (!matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_ = Teuchos::rcp(new ADAPTER::MortarVolCoupl());

    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy =
        Teuchos::rcp(new TSI::UTILS::TSIMaterialStrategy());
    // init coupling adapter projection matrices
    volcoupl_->Init(structdis, thermodis, NULL, NULL, NULL, NULL, materialstrategy);
    // redistribute discretizations to meet needs of volmortar coupling
    volcoupl_->Redistribute();
    // setup projection matrices
    volcoupl_->Setup();
  }

  if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") ==
      INPAR::STR::int_old)
    dserror("old structural time integration no longer supported in tsi");
  else
  {
    Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermo = Teuchos::rcp(
        new ADAPTER::ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams(), thermodis));
    thermo_ = thermo->ThermoFieldrcp();

    //  // access structural dynamic params list which will be possibly modified while creating the
    //  time integrator
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
        ADAPTER::STR::BuildStructureAlgorithm(sdyn);
    adapterbase_ptr->Init(DRT::Problem::Instance()->TSIDynamicParams(),
        const_cast<Teuchos::ParameterList&>(sdyn), structdis);

    // set the temperature; Monolithic does this in it's own constructor with potentially
    // redistributed discretizations
    if (DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(
            DRT::Problem::Instance()->TSIDynamicParams(), "COUPALGO") != INPAR::TSI::Monolithic)
    {
      if (matchinggrid_)
        structdis->SetState(1, "temperature", ThermoField()->Tempnp());
      else
        structdis->SetState(
            1, "temperature", volcoupl_->ApplyVectorMapping12(ThermoField()->Tempnp()));
    }

    adapterbase_ptr->Setup();
    structure_ =
        Teuchos::rcp_dynamic_cast<ADAPTER::StructureWrapper>(adapterbase_ptr->StructureField());

    if (restart &&
        DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(
            DRT::Problem::Instance()->TSIDynamicParams(), "COUPALGO") == INPAR::TSI::Monolithic)
      structure_->Setup();

    StructureField()->Discretization()->ClearState(true);
  }

  // initialise displacement field needed for Output()
  // (get noderowmap of discretisation for creating this multivector)
  // TODO: why nds 0 and not 1????
  dispnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(ThermoField()->Discretization()->NodeRowMap()), 3, true));
  tempnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(StructureField()->Discretization()->NodeRowMap()), 1, true));

  // setup coupling object for matching discretization
  if (matchinggrid_)
  {
    coupST_ = Teuchos::rcp(new ADAPTER::Coupling());
    coupST_->SetupCoupling(*StructureField()->Discretization(), *ThermoField()->Discretization(),
        *StructureField()->Discretization()->NodeRowMap(),
        *ThermoField()->Discretization()->NodeRowMap(), 1, true);
  }

  // setup mortar coupling
  if (DRT::Problem::Instance()->ProblemType() == prb_tsi)
  {
    DRT::Condition* mrtrcond = StructureField()->Discretization()->GetCondition("MortarMulti");
    if (mrtrcond != NULL)
    {
      mortar_coupling_ = Teuchos::rcp(new MORTAR::MultiFieldCoupling());
      mortar_coupling_->PushBackCoupling(
          StructureField()->Discretization(), 0, std::vector<int>(3, 1));
      mortar_coupling_->PushBackCoupling(
          ThermoField()->Discretization(), 0, std::vector<int>(1, 1));
    }
  }

  // reset states
  StructureField()->Discretization()->ClearState(true);
  ThermoField()->Discretization()->ClearState(true);

  return;
}


/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm() {}


/*----------------------------------------------------------------------*
 | update (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField()->Update();
  ThermoField()->Update();
  if (contact_strategy_lagrange_ != Teuchos::null)
    contact_strategy_lagrange_->Update((StructureField()->Dispnp()));

  return;
}


/*----------------------------------------------------------------------*
 | output (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output(bool forced_writerestart)
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("RESULTSEVRY");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");

  //========================
  // output for thermofield:
  //========================
  ApplyStructCouplingState(StructureField()->Dispnp(), StructureField()->Velnp());
  ThermoField()->Output(forced_writerestart);

  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if (forced_writerestart == true and
      ((upres != 0 and (Step() % upres == 0)) or ((uprestart != 0) and (Step() % uprestart == 0))))
  {
    // displacement has already been written into thermo field for this step
    ;
  }
  else if ((upres != 0 and (Step() % upres == 0)) or
           ((uprestart != 0) and (Step() % uprestart == 0)) or forced_writerestart == true)
  {
    if (matchinggrid_)
    {
      OutputDeformationInThr(StructureField()->Dispn(), StructureField()->Discretization());

      ThermoField()->DiscWriter()->WriteVector("displacement", dispnp_, IO::nodevector);
    }
    else
    {
      Teuchos::RCP<const Epetra_Vector> dummy =
          volcoupl_->ApplyVectorMapping21(StructureField()->Dispnp());

      // determine number of space dimensions
      const int numdim = DRT::Problem::Instance()->NDim();

      int err(0);

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (ThermoField()->Discretization()->NumMyRowNodes()); lnodeid++)
      {
        DRT::Node* thermnode = ThermoField()->Discretization()->lRowNode(lnodeid);
        std::vector<int> thermnodedofs_1 = ThermoField()->Discretization()->Dof(1, thermnode);

        // now we transfer displacment dofs only
        for (int index = 0; index < numdim; ++index)
        {
          // global and processor's local fluid dof ID
          const int sgid = thermnodedofs_1[index];
          const int slid = ThermoField()->Discretization()->DofRowMap(1)->LID(sgid);


          // get value of corresponding displacement component
          double disp = (*dummy)[slid];
          // insert velocity value into node-based vector
          err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
          if (err != 0) dserror("error while inserting a value into dispnp_");
        }

        // for security reasons in 1D or 2D problems:
        // set zeros for all unused velocity components
        for (int index = numdim; index < 3; ++index)
        {
          err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
          if (err != 0) dserror("error while inserting a value into dispnp_");
        }
      }  // for lnodid

      ThermoField()->DiscWriter()->WriteVector("displacement", dispnp_, IO::nodevector);
    }
  }


  //===========================
  // output for structurefield:
  //===========================
  ApplyThermoCouplingState(ThermoField()->Tempnp());
  StructureField()->Output(forced_writerestart);

  // mapped temperatures for structure field
  if ((upres != 0 and (Step() % upres == 0)) or ((uprestart != 0) and (Step() % uprestart == 0)) or
      forced_writerestart == true)
    if (not matchinggrid_)
    {
      //************************************************************************************
      Teuchos::RCP<const Epetra_Vector> dummy1 =
          volcoupl_->ApplyVectorMapping12(ThermoField()->Tempnp());

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (StructureField()->Discretization()->NumMyRowNodes());
           lnodeid++)
      {
        DRT::Node* structnode = StructureField()->Discretization()->lRowNode(lnodeid);
        std::vector<int> structdofs = StructureField()->Discretization()->Dof(1, structnode);

        // global and processor's local structure dof ID
        const int sgid = structdofs[0];
        const int slid = StructureField()->Discretization()->DofRowMap(1)->LID(sgid);

        // get value of corresponding displacement component
        double temp = (*dummy1)[slid];
        // insert velocity value into node-based vector
        int err = tempnp_->ReplaceMyValue(lnodeid, 0, temp);
        if (err != 0) dserror("error while inserting a value into tempnp_");
      }  // for lnodid

      StructureField()->Discretization()->Writer()->WriteVector(
          "struct_temperature", tempnp_, IO::nodevector);
    }


  // reset states
  StructureField()->Discretization()->ClearState(true);
  ThermoField()->Discretization()->ClearState(true);
}  // Output()


/*----------------------------------------------------------------------*
 | communicate the displacement vector to THR field          dano 12/11 |
 | enable visualisation of thermal variables on deformed body           |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::OutputDeformationInThr(
    Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<DRT::Discretization> structdis)
{
  if (dispnp == Teuchos::null) dserror("Got null pointer for displacements");

  int err(0);

  // get dofrowmap of structural discretisation
  const Epetra_Map* structdofrowmap = structdis->DofRowMap(0);

  // loop over all local nodes of thermal discretisation
  for (int lnodeid = 0; lnodeid < (ThermoField()->Discretization()->NumMyRowNodes()); lnodeid++)
  {
    // Here we rely on the fact that the thermal discretisation is a clone of
    // the structural mesh.
    // => a thermal node has the same local (and global) ID as its corresponding
    // structural node!

    // get the processor's local structural node with the same lnodeid
    DRT::Node* structlnode = structdis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this structural node
    std::vector<int> structnodedofs = structdis->Dof(0, structlnode);
    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    // now we transfer displacment dofs only
    for (int index = 0; index < numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int sgid = structnodedofs[index];
      const int slid = structdofrowmap->LID(sgid);

      // get value of corresponding displacement component
      double disp = (*dispnp)[slid];
      // insert velocity value into node-based vector
      err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
      if (err != 0) dserror("error while inserting a value into dispnp_");
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index = numdim; index < 3; ++index)
    {
      err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
      if (err != 0) dserror("error while inserting a value into dispnp_");
    }

  }  // for lnodid

  return;

}  // OutputDeformationInThr()


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 12/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> TSI::Algorithm::CalcVelocity(
    Teuchos::RCP<const Epetra_Vector> dispnp)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector(*(StructureField()->Dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *dispnp, -1. / Dt());

  return vel;
}  // CalcVelocity()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ApplyThermoCouplingState(
    Teuchos::RCP<const Epetra_Vector> temp, Teuchos::RCP<const Epetra_Vector> temp_res)
{
  if (matchinggrid_)
  {
    if (temp != Teuchos::null) StructureField()->Discretization()->SetState(1, "temperature", temp);
    if (temp_res != Teuchos::null)
      StructureField()->Discretization()->SetState(1, "residual temperature", temp_res);
  }
  else
  {
    if (temp != Teuchos::null)
      StructureField()->Discretization()->SetState(
          1, "temperature", volcoupl_->ApplyVectorMapping12(temp));
  }
}  // ApplyThermoCouplingState()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ApplyStructCouplingState(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  if (matchinggrid_)
  {
    if (disp != Teuchos::null) ThermoField()->Discretization()->SetState(1, "displacement", disp);
    if (vel != Teuchos::null) ThermoField()->Discretization()->SetState(1, "velocity", vel);
  }
  else
  {
    if (disp != Teuchos::null)
      ThermoField()->Discretization()->SetState(
          1, "displacement", volcoupl_->ApplyVectorMapping21(disp));
    if (vel != Teuchos::null)
      ThermoField()->Discretization()->SetState(
          1, "velocity", volcoupl_->ApplyVectorMapping21(vel));
  }
}  // ApplyStructCouplingState()


/*----------------------------------------------------------------------*/
void TSI::Algorithm::GetContactStrategy()
{
  INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
          DRT::Problem::Instance()->ContactDynamicParams(), "STRATEGY");

  if (stype == INPAR::CONTACT::solution_nitsche)
  {
    if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
            DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") !=
        INPAR::STR::int_standard)
      dserror("thermo-mechanical contact only with new structural time integration");

    if (coupST_ == Teuchos::null) dserror("coupST_ not yet here");

    STR::MODELEVALUATOR::Contact& a = static_cast<STR::MODELEVALUATOR::Contact&>(
        StructureField()->ModelEvaluator(INPAR::STR::model_contact));
    contact_strategy_nitsche_ =
        Teuchos::rcp_dynamic_cast<CONTACT::CoNitscheStrategyTsi>(a.StrategyPtr(), false);
    contact_strategy_nitsche_->EnableRedistribution();

    thermo_->SetNitscheContactStrategy(contact_strategy_nitsche_);

    return;
  }

  else if (stype == INPAR::CONTACT::solution_lagmult)
  {
    if (StructureField()->HaveModel(INPAR::STR::model_contact))
      dserror(
          "structure should not have a Lagrange strategy ... as long as condensed"
          "contact formulations are not moved to the new structural time integration");

    std::vector<DRT::Condition*> ccond(0);
    StructureField()->Discretization()->GetCondition("Contact", ccond);
    if (ccond.size() == 0) return;

    // ---------------------------------------------------------------------
    // create the contact factory
    // ---------------------------------------------------------------------
    CONTACT::STRATEGY::Factory factory;
    factory.Init(Teuchos::rcp_dynamic_cast<DRT::DiscretizationInterface>(
        StructureField()->Discretization(), true));
    factory.Setup();

    // check the problem dimension
    factory.CheckDimension();

    // create some local variables (later to be stored in strategy)
    std::vector<Teuchos::RCP<CONTACT::CoInterface>> interfaces;
    Teuchos::ParameterList cparams;

    // read and check contact input parameters
    factory.ReadAndCheckInput(cparams);

    // ---------------------------------------------------------------------
    // build the contact interfaces
    // ---------------------------------------------------------------------
    // FixMe Would be great, if we get rid of these poro parameters...
    bool poroslave = false;
    bool poromaster = false;
    factory.BuildInterfaces(cparams, interfaces, poroslave, poromaster);

    // ---------------------------------------------------------------------
    // build the solver strategy object
    // ---------------------------------------------------------------------
    contact_strategy_lagrange_ = Teuchos::rcp_dynamic_cast<CONTACT::CoTSILagrangeStrategy>(
        factory.BuildStrategy(cparams, poroslave, poromaster, 1e8, interfaces), true);

    // build the search tree
    factory.BuildSearchTree(interfaces);

    // print final screen output
    factory.Print(interfaces, contact_strategy_lagrange_, cparams);

    // ---------------------------------------------------------------------
    // final touches to the contact strategy
    // ---------------------------------------------------------------------

    contact_strategy_lagrange_->StoreDirichletStatus(StructureField()->GetDBCMapExtractor());

    Teuchos::RCP<Epetra_Vector> zero_disp =
        Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(), true));
    contact_strategy_lagrange_->SetState(MORTAR::state_new_displacement, *zero_disp);
    contact_strategy_lagrange_->SaveReferenceState(zero_disp);
    contact_strategy_lagrange_->EvaluateReferenceState(zero_disp);
    contact_strategy_lagrange_->Inttime_init();
    contact_strategy_lagrange_->RedistributeContact(StructureField()->Dispn());
    contact_strategy_lagrange_->InitBinStrategyforTimestep(StructureField()->Veln());

    if (contact_strategy_lagrange_ != Teuchos::null)
    {
      contact_strategy_lagrange_->SetAlphafThermo(DRT::Problem::Instance()->ThermalDynamicParams());
      contact_strategy_lagrange_->SetCoupling(coupST_);
    }
  }
}
