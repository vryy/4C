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
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"

#include "../drt_adapter/adapter_thermo.H"
#include "../drt_lib/drt_discret.H"

//for coupling of nonmatching meshes
#include "../drt_adapter/adapter_coupling_volmortar.H"
#include "../drt_volmortar/volmortar_utils.H"

// contact
#include "../drt_contact/contact_lagrange_strategy.H"
#include "../drt_contact/contact_tsi_lagrange_strategy.H"
#include "../drt_contact/meshtying_contact_bridge.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(const Epetra_Comm& comm)
: AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
  dispnp_(Teuchos::null),
  tempnp_(Teuchos::null),
  matchinggrid_(DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->TSIDynamicParams(),"MATCHINGGRID")),
  volcoupl_(Teuchos::null)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermodis = DRT::Problem::Instance()->GetDis("thermo");

  if(!matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_=Teuchos::rcp(new ADAPTER::MortarVolCoupl() );

    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy= Teuchos::rcp(new TSI::UTILS::TSIMaterialStrategy() );
    // init coupling adapter projection matrices
    volcoupl_->Init(structdis,thermodis,NULL,NULL,NULL,NULL,materialstrategy);
    // redistribute discretizations to meet needs of volmortar coupling
    volcoupl_->Redistribute();
    // setup projection matrices
    volcoupl_->Setup();
  }

  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure
    = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams(), const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = structure->StructureField();
  structure_->Setup();

  Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermo
    = Teuchos::rcp(new ADAPTER::ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams(),thermodis));
  thermo_ = thermo->ThermoFieldrcp();

  // initialise displacement field needed for Output()
  // (get noderowmap of discretisation for creating this multivector)
  // TODO: why nds 0 and not 1????
  dispnp_ = Teuchos::rcp(new Epetra_MultiVector(*(ThermoField()->Discretization()->NodeRowMap()),3,true));
  tempnp_ = Teuchos::rcp(new Epetra_MultiVector(*(StructureField()->Discretization()->NodeRowMap()),1,true));

  // setup coupling object for matching discretization
  if (matchinggrid_)
  {
    coupST_ = Teuchos::rcp(new ADAPTER::Coupling());
    coupST_->SetupCoupling(*StructureField()->Discretization(),
        *ThermoField()   ->Discretization(),
        *StructureField()->Discretization()->NodeRowMap(),
        *ThermoField()   ->Discretization()->NodeRowMap(),
        1,
        true);
  }

  return;
}


/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*
 | update (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField()->Update();
  ThermoField()->Update();

  // update contact
  if (StructureField()->MeshtyingContactBridge()!=Teuchos::null)
    if (StructureField()->MeshtyingContactBridge()->ContactManager()!=Teuchos::null)
      dynamic_cast<CONTACT::CoTSILagrangeStrategy&>(StructureField()->MeshtyingContactBridge()->
          ContactManager()->GetStrategy()).Update(StructureField()->WriteAccessDispnp(),coupST_);

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

  //===========================
  // output for structurefield:
  //===========================
  ApplyThermoCouplingState(ThermoField()->Tempnp());
  StructureField()->Output(forced_writerestart);

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("RESULTSEVRY");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");

  // mapped temperatures for structure field
  if ( (upres!=0 and (Step()%upres == 0))
    or ( (uprestart != 0) and (Step()%uprestart == 0) )
    or forced_writerestart == true )
    if(not matchinggrid_)
    {
      //************************************************************************************
      Teuchos::RCP<const Epetra_Vector> dummy1 = volcoupl_->ApplyVectorMapping12(ThermoField()->Tempnp());

      // loop over all local nodes of thermal discretisation
      for (int lnodeid=0; lnodeid<(StructureField()->Discretization()->NumMyRowNodes()); lnodeid++)
      {
        DRT::Node* structnode = StructureField()->Discretization()->lRowNode(lnodeid);
        std::vector<int> structdofs = StructureField()->Discretization()->Dof(1,structnode);

        // global and processor's local structure dof ID
        const int sgid = structdofs[0];
        const int slid = StructureField()->Discretization()->DofRowMap(1)->LID(sgid);

        // get value of corresponding displacement component
        double temp = (*dummy1)[slid];
        // insert velocity value into node-based vector
        int err = tempnp_->ReplaceMyValue(lnodeid, 0, temp);
        if (err!= 0) dserror("error while inserting a value into tempnp_");
      } // for lnodid

      StructureField()->DiscWriter()->WriteVector("struct_temperature",tempnp_,IO::DiscretizationWriter::nodevector);
    }

  //========================
  // output for thermofield:
  //========================
  ApplyStructCouplingState(StructureField()->Dispnp(),StructureField()->Velnp());
  ThermoField()->Output(forced_writerestart);

  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if( forced_writerestart == true and
      ( (upres!=0 and (Step()%upres == 0)) or ((uprestart != 0) and (Step()%uprestart == 0)) ) )
  {
    // displacement has already been written into thermo field for this step
    return;
  }
  if ( (upres!=0 and (Step()%upres == 0))
    or ( (uprestart != 0) and (Step()%uprestart == 0) )
    or forced_writerestart == true )
    {
      if(matchinggrid_)
      {
        OutputDeformationInThr(
            StructureField()->Dispn(),
            StructureField()->Discretization()
            );

        ThermoField()->DiscWriter()->WriteVector("displacement",dispnp_,IO::DiscretizationWriter::nodevector);
      }
      else
      {
        Teuchos::RCP<const Epetra_Vector> dummy = volcoupl_->ApplyVectorMapping21(StructureField()->Dispnp());

        // determine number of space dimensions
        const int numdim = DRT::Problem::Instance()->NDim();

        int err(0);

        // loop over all local nodes of thermal discretisation
        for (int lnodeid=0; lnodeid<(ThermoField()->Discretization()->NumMyRowNodes()); lnodeid++)
        {
          DRT::Node* thermnode = ThermoField()->Discretization()->lRowNode(lnodeid);
          std::vector<int> thermnodedofs_1 = ThermoField()->Discretization()->Dof(1,thermnode);

          // now we transfer displacment dofs only
          for(int index=0; index<numdim; ++index)
          {
            // global and processor's local fluid dof ID
            const int sgid = thermnodedofs_1[index];
            const int slid = ThermoField()->Discretization()->DofRowMap(1)->LID(sgid);


            // get value of corresponding displacement component
            double disp = (*dummy)[slid];
            // insert velocity value into node-based vector
            err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
            if (err!= 0) dserror("error while inserting a value into dispnp_");
          }

          // for security reasons in 1D or 2D problems:
          // set zeros for all unused velocity components
          for (int index=numdim; index < 3; ++index)
          {
            err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
            if (err!= 0) dserror("error while inserting a value into dispnp_");
          }
        } // for lnodid

        ThermoField()->DiscWriter()->WriteVector("displacement",dispnp_,IO::DiscretizationWriter::nodevector);
      }
    }

  //reset states
  StructureField()->Discretization()->ClearState(true);
  ThermoField()->Discretization()->ClearState(true);
}  // Output()


/*----------------------------------------------------------------------*
 | communicate the displacement vector to THR field          dano 12/11 |
 | enable visualisation of thermal variables on deformed body           |
 *----------------------------------------------------------------------*/
void  TSI::Algorithm::OutputDeformationInThr(
  Teuchos::RCP<const Epetra_Vector> dispnp,
  Teuchos::RCP<DRT::Discretization> structdis
  )
{
  if (dispnp == Teuchos::null)
    dserror("Got null pointer for displacements");

  int err(0);

  // get dofrowmap of structural discretisation
  const Epetra_Map* structdofrowmap = structdis->DofRowMap(0);

  // loop over all local nodes of thermal discretisation
  for (int lnodeid=0; lnodeid<(ThermoField()->Discretization()->NumMyRowNodes()); lnodeid++)
  {
    // Here we rely on the fact that the thermal discretisation is a clone of
    // the structural mesh.
    // => a thermal node has the same local (and global) ID as its corresponding
    // structural node!

    // get the processor's local structural node with the same lnodeid
    DRT::Node* structlnode = structdis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this structural node
    std::vector<int> structnodedofs = structdis->Dof(0,structlnode);
    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    // now we transfer displacment dofs only
    for(int index=0; index<numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int sgid = structnodedofs[index];
      const int slid = structdofrowmap->LID(sgid);

      // get value of corresponding displacement component
      double disp = (*dispnp)[slid];
      // insert velocity value into node-based vector
      err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
      if (err!= 0) dserror("error while inserting a value into dispnp_");
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index=numdim; index < 3; ++index)
    {
      err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
      if (err!= 0) dserror("error while inserting a value into dispnp_");
    }

  } // for lnodid

  return;

}  // OutputDeformationInThr()


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 12/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> TSI::Algorithm::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> dispnp
  )
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector( *(StructureField()->Dispn()) ) );
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1./Dt(), *dispnp, -1./Dt());

  return vel;
}  // CalcVelocity()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ApplyThermoCouplingState(Teuchos::RCP<const Epetra_Vector> temp,
                                              Teuchos::RCP<const Epetra_Vector> temp_res)
{
  if (matchinggrid_)
  {
    if (temp != Teuchos::null)
      StructureField()->Discretization()->SetState(1,"temperature",temp);
    if (temp_res != Teuchos::null)
      StructureField()->Discretization()->SetState(1,"residual temperature",temp_res);
  }
  else
  {
    if (temp != Teuchos::null)
      StructureField()->Discretization()->SetState(1,"temperature",volcoupl_->ApplyVectorMapping12(temp));
  }
}  // ApplyThermoCouplingState()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ApplyStructCouplingState(Teuchos::RCP<const Epetra_Vector> disp,
                                              Teuchos::RCP<const Epetra_Vector> vel)
{
  if (matchinggrid_)
  {
    if (disp != Teuchos::null)
      ThermoField()->Discretization()->SetState(1, "displacement", disp);
    if (vel != Teuchos::null)
      ThermoField()->Discretization()->SetState(1, "velocity", vel);
  }
  else
  {
    if (disp != Teuchos::null)
      ThermoField()->Discretization()->SetState(1,"displacement",volcoupl_->ApplyVectorMapping21(disp));
    if (vel != Teuchos::null)
      ThermoField()->Discretization()->SetState(1,"velocity",volcoupl_->ApplyVectorMapping21(vel));
  }
}  // ApplyStructCouplingState()


/*----------------------------------------------------------------------*/
