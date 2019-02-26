/*----------------------------------------------------------------------*/
/*!
\file xfluid_state_creator.cpp

\brief Creates a state object for (in)stationary XFEM fluid problems

\level 0

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/


#include "xfluid_state_creator.H"

#include "xfluid_state.H"
#include "xfluidfluid_state.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../linalg/linalg_mapextractor.H"

#include "../drt_cut/cut_cutwizard.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_dofset.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 |  Perform the cut and fill state container               schott 01/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidState> FLD::XFluidStateCreator::Create(
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,  //!< xfluid background discretization
    Teuchos::RCP<const Epetra_Vector>
        back_disp_col,  //!< col vector holding background ALE displacements for backdis
    Teuchos::ParameterList& solver_params,  //!< solver parameters
    const int step,                         //!< current time step
    const double& time                      //!< current time
)
{
#ifdef DEBUG
  if (condition_manager_ == Teuchos::null) dserror("no condition manager available!");
#endif

  //--------------------------------------------------------------------------------------
  // create new cut wizard &dofset
  Teuchos::RCP<GEO::CutWizard> wizard;
  Teuchos::RCP<XFEM::XFEMDofSet> dofset;

  CreateNewCutState(dofset, wizard, xdiscret, back_disp_col, solver_params, step);

  //--------------------------------------------------------------------------------------
  // Create the XFluidState object

  Teuchos::RCP<const Epetra_Map> xfluiddofrowmap =
      Teuchos::rcp(new Epetra_Map(*xdiscret->DofRowMap()));

  Teuchos::RCP<const Epetra_Map> xfluiddofcolmap =
      Teuchos::rcp(new Epetra_Map(*xdiscret->DofColMap()));

  Teuchos::RCP<XFluidState> state = Teuchos::rcp(
      new FLD::XFluidState(condition_manager_, wizard, dofset, xfluiddofrowmap, xfluiddofcolmap));

  //--------------------------------------------------------------------------------------
  state->SetupMapExtractors(xdiscret, time);

  return state;
}


/*----------------------------------------------------------------------*
 |  Perform the cut and fill state container                kruse 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidFluidState> FLD::XFluidStateCreator::Create(
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,     //!< xfluid background discretization
    const Teuchos::RCP<DRT::Discretization>& embfluiddiscret,  //!< embedded fluid discretization
    Teuchos::RCP<const Epetra_Vector>
        back_disp_col,  //!< col vector holding background ALE displacements for backdis
    Teuchos::ParameterList& solver_params,  //!< solver parameters
    const int step,                         //!< current time step
    const double& time                      //!< current time
)
{
#ifdef DEBUG
  if (condition_manager_ == Teuchos::null) dserror("no condition manager available!");
#endif

  //--------------------------------------------------------------------------------------
  // create new cut wizard & dofset
  Teuchos::RCP<GEO::CutWizard> wizard;
  Teuchos::RCP<XFEM::XFEMDofSet> dofset;

  CreateNewCutState(dofset, wizard, xdiscret, back_disp_col, solver_params, step);

  //--------------------------------------------------------------------------------------
  // Create the XFluidFluidState object

  Teuchos::RCP<const Epetra_Map> xfluiddofrowmap =
      Teuchos::rcp(new Epetra_Map(*xdiscret->DofRowMap()));

  Teuchos::RCP<const Epetra_Map> xfluiddofcolmap =
      Teuchos::rcp(new Epetra_Map(*xdiscret->DofColMap()));

  Teuchos::RCP<const Epetra_Map> embfluiddofrowmap =
      Teuchos::rcp(new Epetra_Map(*embfluiddiscret->DofRowMap()));

  Teuchos::RCP<FLD::XFluidFluidState> state = Teuchos::rcp(new FLD::XFluidFluidState(
      condition_manager_, wizard, dofset, xfluiddofrowmap, xfluiddofcolmap, embfluiddofrowmap));

  //--------------------------------------------------------------------------------------
  state->SetupMapExtractors(xdiscret, embfluiddiscret, time);

  return state;
}


/*----------------------------------------------------------------------*
 |  Initialize ALE state vectors                           schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::CreateNewCutState(
    Teuchos::RCP<XFEM::XFEMDofSet>& dofset,  //!< xfem dofset obtained from the new wizard
    Teuchos::RCP<GEO::CutWizard>&
        wizard,  //!< cut wizard associated with current intersection state
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,  //!< xfluid background discretization
    Teuchos::RCP<const Epetra_Vector>
        back_disp_col,  //!< col vector holding background ALE displacements for backdis
    Teuchos::ParameterList& solver_params,  //!< solver parameters
    const int step                          //!< current time step
)
{
  // new wizard using information about cutting sides from the condition_manager
  wizard = Teuchos::rcp(new GEO::CutWizard(xdiscret));

  // Set options for the cut wizard
  wizard->SetOptions(nodal_dofset_strategy_,  // strategy for nodal dofset management
      VolumeCellGaussPointBy_,                // how to create volume cell Gauss points?
      BoundCellGaussPointBy_,                 // how to create boundary cell Gauss points?
      gmsh_cut_out_,                          // gmsh output for cut library
      true,                                   // find point positions
      false,                                  // generate only tet cells
      true                                    // print screen output
  );

  //--------------------------------------------------------------------------------------
  // set state for all mesh cutting

  // loop all mesh coupling objects
  for (int mc_idx = 0; mc_idx < condition_manager_->NumMeshCoupling(); mc_idx++)
  {
    Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx);

    if (!mc_coupl->CutGeometry()) continue;  // If don't cut the background mesh.

    wizard->AddCutterState(mc_idx, mc_coupl->GetCutterDis(), mc_coupl->GetCutterDispCol(),
        condition_manager_->GetMeshCouplingStartGID(mc_idx));
  }

  //--------------------------------------------------------------------------------------
  // set background state (background mesh displacements and level-set values)

  wizard->SetBackgroundState(
      back_disp_col,  //!< col vector holding background ALE displacements for backdis
      condition_manager_
          ->GetLevelSetFieldCol(),  //!< col vector holding nodal level-set values based on backdis
      condition_manager_->GetLevelSetCouplingGid()  //!< global side id for level-set coupling
  );

  //--------------------------------------------------------------------------------------
  // Initialize cut objects into the cut
  wizard->Prepare();

  // Loop all mesh coupling objects:
  // -- Find corresponding marked surfaces loaded into the cut.
  for (int mc_idx = 0; mc_idx < condition_manager_->NumMeshCoupling(); mc_idx++)
  {
    Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx);

    if (mc_coupl->IsMarkedGeometry())
    {
      wizard->SetMarkedConditionSides(
          mc_coupl->GetCutterDis(), condition_manager_->GetMeshCouplingStartGID(mc_idx));
    }
  }

  //--------------------------------------------------------------------------------------
  // performs the "CUT"
  wizard->Cut(include_inner_);

  //--------------------------------------------------------------------------------------
  // set the new dofset after cut
  int maxNumMyReservedDofsperNode = (maxnumdofsets_)*4;

  // create a new XFEM-dofset
  dofset = Teuchos::rcp(new XFEM::XFEMDofSet(*wizard, maxNumMyReservedDofsperNode, *xdiscret));

  const int restart = DRT::Problem::Instance()->Restart();
  if ((step < 1) or restart) minnumdofsets_ = xdiscret->DofRowMap()->MinAllGID();

  dofset->SetMinGID(minnumdofsets_);         // set the minimal GID of xfem dis
  xdiscret->ReplaceDofSet(0, dofset, true);  // fluid dofset has nds = 0

  xdiscret->FillComplete(true, false, false);

  // print all dofsets
  xdiscret->GetDofSetProxy()->PrintAllDofsets(xdiscret->Comm());

  //--------------------------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  // REMARK: this has to be done after replacing the discret' dofset (via discret_->ReplaceDofSet)
  xdiscret->ComputeNullSpaceIfNecessary(solver_params, true);
}
