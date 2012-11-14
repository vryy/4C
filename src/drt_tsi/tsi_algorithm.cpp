/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.cpp

\brief  Basis of all TSI algorithms that perform a coupling between the linear
        momentum equation and the heat conduction equation

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
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
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(const Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams())
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams(), structdis));
  structure_ = structure->StructureFieldrcp();

  Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermo = Teuchos::rcp(new ADAPTER::ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()));
  thermo_ = thermo->ThermoFieldrcp();

  // initialise displacement field needed for Output()
  // (get noderowmap of discretisation for creating this multivector)
  dispnp_ = Teuchos::rcp(new Epetra_MultiVector(*(ThermoField()->Discretization()->NodeRowMap()),3,true));

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
  return;
}


/*----------------------------------------------------------------------*
 | output (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output()
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.
  StructureField()->Output();

  ThermoField()->Output();

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("UPRES");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");
  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if ( (upres!=0 and (Step()%upres == 0))
    or ( (uprestart != 0) and (Step()%uprestart == 0) ) )
    {
      OutputDeformationInThr(
          StructureField()->Dispn(),
          StructureField()->Discretization()
          );

      ThermoField()->DiscWriter()->WriteVector("displacement",dispnp_,IO::DiscretizationWriter::nodevector);
    }

}  // Output


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
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> dispnp
  )
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector( *(StructureField()->ExtractDispn()) ) );
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1./Dt(), *dispnp, -1./Dt());

  return vel;
}  // CalcVelocity()


/*----------------------------------------------------------------------*/
