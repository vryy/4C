/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_utils.cpp

 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "poromultiphase_utils.H"

#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "poromultiphase_utils_clonestrategy.H"

#include "poromultiphase_partitioned.H"
#include "poromultiphase_partitioned_twoway.H"

#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele.H"
#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_lib/drt_utils_createdis.H"

/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    int& nds_disp,
    int& nds_vel,
    int& nds_solidpressure)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);
  if(!structdis->Filled())
    structdis->FillComplete();
  if(!fluiddis->Filled())
    fluiddis->FillComplete();

  if (fluiddis->NumGlobalNodes()==0)
  {
    // fill poro fluid discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy>(
        structdis,
        fluiddis);
  }
  else
  {
    dserror("Fluid discretization given in input file. This is not supported!");
  }

  structdis->FillComplete();
  fluiddis->FillComplete();

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();

  //assign structure dof set to fluid and save the dofset number
  nds_disp = fluiddis->AddDofSet(structdofset);
  if (nds_disp!=1)
    dserror("unexpected dof sets in porofluid field");
  // velocities live on same dofs as displacements
  nds_vel=nds_disp;

  if (structdis->AddDofSet(fluiddofset)!=1)
    dserror("unexpected dof sets in structure field");

  // build auxiliary dofset for postprocessing solid pressures
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(1,0,0,false));
  nds_solidpressure = fluiddis->AddDofSet(dofsetaux);
  // add it also to the solid field
  structdis->AddDofSet(fluiddis->GetDofSetProxy(nds_solidpressure));

  structdis->FillComplete();
  fluiddis->FillComplete();

  return;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::AssignMaterialPointers(
    const std::string& struct_disname,
    const std::string& fluid_disname)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);

  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,fluiddis);
}

/*----------------------------------------------------------------------*
 | setup algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASE::PoroMultiPhaseBase> POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(
    INPAR::POROMULTIPHASE::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<POROMULTIPHASE::PoroMultiPhaseBase> algo = Teuchos::null;

  switch(solscheme)
  {
  case INPAR::POROMULTIPHASE::solscheme_twoway:
  {
    // call constructor
    algo =
        Teuchos::rcp(new POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay(
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

