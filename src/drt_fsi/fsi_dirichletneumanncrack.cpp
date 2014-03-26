/*
 * fsidirichletneumanncrack.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: sudhakar
 */

#include "fsi_dirichletneumanncrack.H"

#include "../drt_adapter/ad_str_fsi_crack.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_xfem.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*-------------------------------------------------------------------------------------*
 *        Constructor
 *-------------------------------------------------------------------------------------*/
FSI::DirichletNeumann_Crack::DirichletNeumann_Crack(const Epetra_Comm& comm)
  : DirichletNeumann(comm)
{
  crackUpdate_ = false;
}


/*-------------------------------------------------------------------------------------*
 * Central routine that calls all the necessary operations : build interface
 * from fluid side, from structure side, and perform map extraction          sudhakar 03/14
 *-------------------------------------------------------------------------------------*/
void FSI::DirichletNeumann_Crack::update_FSI_interface_Crack()
{
  // propagate crack within the structure
  bool isCrackProp = StructureField()->UpdateCrackInformation( StructureField()->Dispnp() );

  if( isCrackProp )
  {
    // add newly created cut surfaces to the XFEM cut interface
    AddNewCrackSurfaceToCutInterface();

    // rebuild FSI interface of fluid side
    MBFluidField().RebuildFSIInterface();

    // rebuild FSI interface of structure side
    StructureField()->RebuildInterface();

    // Setup coupling conditions on the updated interface
    ADAPTER::Coupling& coupsf = StructureFluidCoupling();
    const int ndim = DRT::Problem::Instance()->NDim();
    coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                   StructureField()->Interface()->FSICondMap(),
                                   *MBFluidField().Discretization(),
                                   MBFluidField().Interface().FSICondMap(),
                                  "FSICoupling",
                                  ndim);

    if (coupsf.MasterDofMap()->NumGlobalElements()==0)
      dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");
  }
}

/*-------------------------------------------------------------------------------------*
 * When crack propagates, it creates new crack surfaces. These new surfaces
 * are added to the XFEM cut discretization. Also, the EpetraVectors that       sudhakar 03/14
 * are defined on the boundary discretization are rebuild to account for
 * the added new node
 *-------------------------------------------------------------------------------------*/
void FSI::DirichletNeumann_Crack::AddNewCrackSurfaceToCutInterface()
{
  const Teuchos::RCP<ADAPTER::FSICrackingStructure>& structfield =
                                Teuchos::rcp_dynamic_cast<ADAPTER::FSICrackingStructure>(StructureField());

  ADAPTER::FluidXFEM& ad_xfem = dynamic_cast<ADAPTER::FluidXFEM&>(MBFluidField());
  ADAPTER::Fluid& ad_flui = ad_xfem.FluidField();

  Teuchos::RCP<DRT::Discretization> boundary_dis = Teuchos::null;
  ad_flui.BoundaryDis( boundary_dis );

  std::map<int, LINALG::Matrix<3,1> > tip_nodes;
  ad_flui.GetCrackTipNodes( tip_nodes );
  structfield->addCrackSurfacesToCutSides( boundary_dis, tip_nodes );
  ad_flui.setBoundaryDis( boundary_dis );
  ad_flui.SetCrackTipNodes( tip_nodes );

  ad_flui.UpdateBoundaryValuesAfterCrack( structfield->getOldNewCrackNodes() );
}
