/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_utils.cpp

 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_utils.H"

#include "poromultiphase_scatra_partitioned_twoway.H"

#include "../drt_poromultiphase/poromultiphase_utils.H"

#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../drt_inpar/inpar_poromultiphase_scatra.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_utils_createdis.H"

/*----------------------------------------------------------------------*
 | setup algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase>
POROMULTIPHASESCATRA::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
    INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase> algo;

  switch(solscheme)
  {
  case INPAR::POROMULTIPHASESCATRA::solscheme_twoway:
  {
    // call constructor
    algo =
        Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay(
                comm,
                timeparams));
    break;
  }
  default:
    dserror("Unknown time-integration scheme for multiphase poro fluid problem");
    break;
  }

  return algo;
}


/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::UTILS::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    int& ndsporo_disp,
    int& ndsporo_vel,
    int& ndsporo_solidpressure,
    int& ndsporofluid_scatra)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.
  //            Then, a scatra disc. is cloned.

  POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
      comm,
      struct_disname,
      fluid_disname,
      ndsporo_disp,
      ndsporo_vel,
      ndsporo_solidpressure);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  // fill scatra discretization by cloning structure discretization
  DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,scatradis);
  scatradis->FillComplete();
  // set implementation type
  for(int i=0; i<scatradis->NumMyColElements(); ++i)
  {
    DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
    if(element == NULL)
      dserror("Invalid element type!");
    else
      element->SetImplType(DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(problem->PoroMultiPhaseScatraDynamicParams(),"SCATRATYPE"));
  }

  // the problem is two way coupled, thus each discretization must know the other discretization

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure/fluid field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // check if ScatraField has 2 discretizations, so that coupling is possible
  if (scatradis->AddDofSet(structdofset) != 1)
    dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddofset) != 2)
    dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddis->GetDofSetProxy(ndsporo_solidpressure)) != 3)
    dserror("unexpected dof sets in scatra field");
  if (structdis->AddDofSet(scatradofset)!=3)
    dserror("unexpected dof sets in structure field");

  ndsporofluid_scatra = fluiddis->AddDofSet(scatradofset);
  if (ndsporofluid_scatra!=3)
    dserror("unexpected dof sets in fluid field");

  structdis->FillComplete(true,false,false);
  fluiddis->FillComplete(true,false,false);
  scatradis->FillComplete(true,false,false);

  return;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::UTILS::AssignMaterialPointers(
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname)
{
  POROMULTIPHASE::UTILS::AssignMaterialPointers(
      struct_disname,
      fluid_disname);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,scatradis);
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(fluiddis,scatradis);
}
