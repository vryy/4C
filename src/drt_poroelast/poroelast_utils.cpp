/*----------------------------------------------------------------------*/
/*!
 \file tsi_utils.cpp

 \brief utility functions for poroelasticity problems

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
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"
#include <Epetra_Time.h>

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
    dserror("Only fluid-poro material type allowed for poroelasticity. Cannot generate fluid mesh.");

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
  // we only support fluid elements here
  eletype.push_back("FLUID3");

  return true; // yes, we copy EVERY element (no submeshes)
}

/*----------------------------------------------------------------------*
 | setup Poroelasticity                                                            |
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
  if (fluiddis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);
  }
  else
    dserror("Structure AND Fluid discretization present. This is not supported.");
  }

  /*----------------------------------------------------------------------*/
