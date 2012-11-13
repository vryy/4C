/*----------------------------------------------------------------------*/
/*!
 \file tsi_utils.cpp

 \brief utility functions for porous media problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 */

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "poroelast_utils.H"
#include "poro_base.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Epetra_Time.h>
#include <Epetra_MpiComm.h>

#include "poro_partitioned.H"
#include "poroelast_monolithic.H"
#include "poro_monolithicstructuresplit.H"
#include "poro_monolithicfluidsplit.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_so3/so3_poro.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(pair<std::string, std::string> ("PoroDirichlet",
      "Dirichlet"));
  conditions_to_copy.insert(pair<std::string, std::string> ("PoroPointNeumann",
      "PointNeumann"));
  conditions_to_copy.insert(pair<std::string, std::string> ("PoroLineNeumann",
      "LineNeumann"));
  conditions_to_copy.insert(pair<std::string, std::string> (
      "PoroSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(pair<std::string, std::string> (
      "PoroVolumeNeumann", "VolumeNeumann"));

  conditions_to_copy.insert(pair<std::string, std::string> ("NoPenetration",
      "NoPenetration"));
  conditions_to_copy.insert(pair<std::string, std::string> ("PoroPartInt",
      "PoroPartInt"));
  conditions_to_copy.insert(pair<std::string, std::string> ("PoroCoupling",
      "PoroCoupling"));
  conditions_to_copy.insert(pair<std::string, std::string> ("FSICoupling",
      "FSICoupling"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(
    const int matid)
{
  //  //// We take the material with the ID specified by the user
  //  //// Here we check first, whether this material is of admissible type
  //  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  //  if ((mtype != INPAR::MAT::m_fluidporo))
  //  dserror("Material with ID %d is not admissible for fluid poroelasticity elements",matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(Teuchos::RCP<
    DRT::Element> newele, DRT::Element* oldele, const int matid,
    const bool isnurbs)
{
  // We must not add a new material type here because that might move
  // the internal material vector. And each element material might
  // have a pointer to that vector. Too bad.
  // So we search for a FLuid material and take the first one we find.
  // => matid from outside remains unused!
  const int matnr = DRT::Problem::Instance()->Materials()->FirstIdByType(
      INPAR::MAT::m_fluidporo);
  if (matnr == -1)
    dserror("Only fluid-poro material type allowed for deformable porous media. Cannot generate fluid mesh.");

    // We need to set material and possibly other things to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property

    //RCP<MAT::Material > mat = oldele->Material();
    //const int matnr = (mat->Parameter()->Id())+1;

#ifdef D_FLUID3 
        DRT::ELEMENTS::Fluid* fluid = dynamic_cast<DRT::ELEMENTS::Fluid*>(newele.get());
        if (fluid!=NULL)
        {
          fluid->SetMaterial(matnr);
          fluid->SetDisType(oldele->Shape()); // set distype as well!
          fluid->SetIsAle(true);
        }
        else
#endif
        {
          dserror("unsupported element type '%s'", typeid(*newele).name());
        }
        return;
      }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, vector<string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support fluid elements here
    eletype.push_back("FLUID3");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastCloneStrategy::CheckPoro(
    DRT::Element* actele)
{
  //all poro elements need to be listed here

  //check for hex8
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* poroelehex8 =
      dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(actele);
  if (poroelehex8!=NULL)
    return true;

  //check for tet4
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* poroeletet4 =
      dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(actele);
  if (poroeletet4!=NULL)
    return true;

  //check for hex27
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* poroelehex27 =
      dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(actele);
  if (poroelehex27!=NULL)
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 | setup Poro discretization                                            |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::SetupPoro(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = problem->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
    structdis->FillComplete();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
  fluiddis = problem->GetDis("fluid");
  if (!fluiddis->Filled())
    fluiddis->FillComplete();

  // we use the structure discretization as layout for the fluid discretization
  if (structdis->NumGlobalNodes() == 0)
    dserror("Structure discretization is empty!");

  // create fluid elements if the fluid discretization is empty
  // create fluid elements if the fluid discretization is empty
  if (fluiddis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);
  }
  else
    dserror("Structure AND Fluid discretization present. This is not supported.");
}


/*----------------------------------------------------------------------*
 | setup Poro algorithm                                            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROELAST::PoroBase> POROELAST::UTILS::CreatePoroAlgorithm(
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams =
      problem->StructuralDynamicParams();

  std::string damping = sdynparams.get<std::string>("DAMPING");
  if(damping != "Material")
    dserror("Material damping has to be used for porous media! Set DAMPING to Material in the STRUCTURAL DYNAMIC section.");

  const Teuchos::ParameterList& poroelastdyn  = problem->PoroelastDynamicParams();

  const INPAR::POROELAST::SolutionSchemeOverFields coupling =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
          poroelastdyn, "COUPALGO");

  // create an empty Poroelast::Algorithm instance
  Teuchos::RCP<POROELAST::PoroBase> poroalgo = Teuchos::null;

  switch (coupling)
  {
    case INPAR::POROELAST::Monolithic:
    {
      // create an POROELAST::Monolithic instance
      poroalgo = Teuchos::rcp(new POROELAST::Monolithic(comm, timeparams));
      break;
    } // monolithic case
    case INPAR::POROELAST::Monolithic_structuresplit:
    {
      // create an POROELAST::MonolithicStructureSplit instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicStructureSplit(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Monolithic_fluidsplit:
    {
      // create an POROELAST::MonolithicFluidSplit instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicFluidSplit(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Partitioned:
    {
      // create an POROELAST::Monolithic instance
      poroalgo = Teuchos::rcp(new POROELAST::Partitioned(comm, timeparams));
      break;
    } // monolithic case
    default:
      dserror("Unknown solutiontype for poroelasticity: %d",coupling);
      break;
  } // end switch

  return poroalgo;
}
