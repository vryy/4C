/*----------------------------------------------------------------------*/
/*!
\file fsi_xfem_algorithm.cpp

\brief Basis of monolithic XFSI algorithm that performs a coupling between the
       structural field equation and XFEM fluid field equations

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                             schott 08/13 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                 schott 08/13 |
 *----------------------------------------------------------------------*/
#include "fsi_xfem_algorithm.H"

//#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_inpar/inpar_fsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                    schott 08/13 |
 *----------------------------------------------------------------------*/
FSI::AlgorithmXFEM::AlgorithmXFEM(const Epetra_Comm& comm,
                                  const Teuchos::ParameterList& timeparams)
  : AlgorithmBase(comm, timeparams)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());


  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,false));
  fluid_ = fluid->FluidFieldrcp();





  //TODO:?
//  // initialise displacement field needed for Output()
//  // (get noderowmap of discretisation for creating this multivector)
//  dispnp_ = Teuchos::rcp(new Epetra_MultiVector(*(ThermoField()->Discretization()->NodeRowMap()),3,true));

  return;
}


/*----------------------------------------------------------------------*
 | destructor (public)                                     schott 08/13 |
 *----------------------------------------------------------------------*/
FSI::AlgorithmXFEM::~AlgorithmXFEM()
{
}


/*----------------------------------------------------------------------*
 | update (protected)                                      schott 08/13 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::Update()
{
  StructureField()->Update();
  FluidField()->Update();
  return;
}


/*----------------------------------------------------------------------*
 | output (protected)                                      schott 08/13 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::Output()
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.
  StructureField()->Output();

  FluidField()->Output();

  // call the FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = fsidyn.get<int>("UPRES");
  int uprestart = fsidyn.get<int>("RESTARTEVRY");
  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if ( (upres!=0 and (Step()%upres == 0))
    or ( (uprestart != 0) and (Step()%uprestart == 0) ) )
    {
      OutputDeformationInThr(
          StructureField()->Dispn(),
          StructureField()->Discretization()
          );

      FluidField()->DiscWriter()->WriteVector("displacement",dispnp_,IO::DiscretizationWriter::nodevector);
    }

}  // Output()


/*----------------------------------------------------------------------*
 | communicate the displacement vector to THR field          dano 12/11 |
 | enable visualisation of thermal variables on deformed body           |
 *----------------------------------------------------------------------*/
void  FSI::AlgorithmXFEM::OutputDeformationInThr(
  Teuchos::RCP<const Epetra_Vector> dispnp,
  Teuchos::RCP<DRT::Discretization> structdis
  )
{
//  if (dispnp == Teuchos::null)
//    dserror("Got null pointer for displacements");
//
//  int err(0);
//
//  // get dofrowmap of structural discretisation
//  const Epetra_Map* structdofrowmap = structdis->DofRowMap(0);
//
//  // loop over all local nodes of thermal discretisation
//  for (int lnodeid=0; lnodeid<(ThermoField()->Discretization()->NumMyRowNodes()); lnodeid++)
//  {
//    // Here we rely on the fact that the thermal discretisation is a clone of
//    // the structural mesh.
//    // => a thermal node has the same local (and global) ID as its corresponding
//    // structural node!
//
//    // get the processor's local structural node with the same lnodeid
//    DRT::Node* structlnode = structdis->lRowNode(lnodeid);
//    // get the degrees of freedom associated with this structural node
//    std::vector<int> structnodedofs = structdis->Dof(0,structlnode);
//    // determine number of space dimensions
//    const int numdim = DRT::Problem::Instance()->NDim();
//
//    // now we transfer displacment dofs only
//    for(int index=0; index<numdim; ++index)
//    {
//      // global and processor's local fluid dof ID
//      const int sgid = structnodedofs[index];
//      const int slid = structdofrowmap->LID(sgid);
//
//      // get value of corresponding displacement component
//      double disp = (*dispnp)[slid];
//      // insert velocity value into node-based vector
//      err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
//      if (err!= 0) dserror("error while inserting a value into dispnp_");
//    }
//
//    // for security reasons in 1D or 2D problems:
//    // set zeros for all unused velocity components
//    for (int index=numdim; index < 3; ++index)
//    {
//      err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
//      if (err!= 0) dserror("error while inserting a value into dispnp_");
//    }
//
//  } // for lnodid

return;

}  // OutputDeformationInThr()


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 12/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::AlgorithmXFEM::CalcVelocity(
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
