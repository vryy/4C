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
#include "poromultiphase_utils_clonestrategy.H"

#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele.H"


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

    // set implementation type
    for(int i=0; i<fluiddis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::PoroFluidMultiPhase* element = dynamic_cast<DRT::ELEMENTS::PoroFluidMultiPhase*>(fluiddis->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");
    }
  }
  else
  {
    dserror("Fluid discretization given in input file. This is not supported!");
  }

  structdis->FillComplete();
  fluiddis->FillComplete();

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> fluiddofset = fluiddis->GetDofSetProxy();

  //assign structure dof set to fluid and save the dofset number
  nds_disp = fluiddis->AddDofSet(structdofset);
  if (nds_disp!=1)
    dserror("unexpected dof sets in porofluid field");
  // velocities live on same dofs as displacements
  nds_vel=nds_disp;

  if (structdis->AddDofSet(fluiddofset)!=1)
    dserror("unexpected dof sets in structure field");

  // build auxiliary dofset for postprocessing solid pressures
  nds_solidpressure = fluiddis->BuildDofSetAuxProxy(1,0,0,false);
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

  const int numelements = fluiddis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* fluidtele = fluiddis->lColElement(i);
    const int gid = fluidtele->Id();

    DRT::Element* structele = structdis->gElement(gid);

    //for coupling we add the source material to the target element and vice versa
    fluidtele->AddMaterial(structele->Material());
    structele->AddMaterial(fluidtele->Material());
  }
}

