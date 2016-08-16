/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumanncrack.cpp

\brief Dirichlet-Neumann partitioning XFSI approach to simulate FSI
       with cracking structures

\maintainer Sudhakar Y.

\level 2
*/
/*----------------------------------------------------------------------*/

#include "fsi_dirichletneumanncrack.H"

#include "../drt_adapter/ad_str_fsi_crack.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_xfsi.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_fluid_xfem.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_xfluid/xfluid.H"

#include "../drt_xfem/xfem_condition_manager.H"

/*-------------------------------------------------------------------------------------*
 *        Constructor
 *-------------------------------------------------------------------------------------*/
FSI::DirichletNeumann_Crack::DirichletNeumann_Crack(const Epetra_Comm& comm)
: DirichletNeumann(comm),
  crackUpdate_(false)
{
  // empty constructor
}


void FSI::DirichletNeumann_Crack::Setup()
{
  // call setup of base class
  FSI::DirichletNeumann::Setup();
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
    Teuchos::RCP<ADAPTER::XFluidFSI> xfluidfsi = Teuchos::rcp_dynamic_cast<ADAPTER::XFluidFSI>(MBFluidField()->FluidField());

    xfluidfsi->RebuildFSIStructInterface();

    // rebuild FSI interface of structure side
    StructureField()->RebuildInterface();

    // Setup coupling conditions on the updated interface
    ADAPTER::Coupling& coupsf = StructureFluidCoupling();
    const int ndim = DRT::Problem::Instance()->NDim();
    coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                   StructureField()->Interface()->FSICondMap(),
                                   *xfluidfsi->BoundaryDiscretization(),
                                   xfluidfsi->StructInterface()->FSICondMap(),
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

  const std::string condname = "XFEMSurfCrackFSIPart";

  Teuchos::RCP<FLD::XFluid> xfluid =  Teuchos::rcp_dynamic_cast<ADAPTER::XFluidFSI>(MBFluidField()->FluidField())->MyFluid();
  Teuchos::RCP<XFEM::MeshCouplingFSICrack> coupl = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSICrack>(xfluid->GetMeshCoupling(condname));

  Teuchos::RCP<DRT::Discretization> boundary_dis = coupl->GetCutterDis();

  std::map<int, LINALG::Matrix<3,1> > & tip_nodes = coupl->GetCrackTipNodes();

  // create new boundary discretization if necessary
  structfield->addCrackSurfacesToCutSides( boundary_dis, tip_nodes );
  if(boundary_dis == Teuchos::null)
    dserror( "Boundary discretization can't be empty" );

  coupl->SetCutterDis( boundary_dis );
  coupl->SetCrackTipNodes( tip_nodes );
  coupl->UpdateBoundaryValuesAfterCrack( structfield->getOldNewCrackNodes() );
}
