/*----------------------------------------------------------------------*/
/*!
 \file poro_utils_clonestrategy.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poroelast_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_fluid_ele/fluid_ele_poro.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "poro_utils_clonestrategy.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroDirichlet",
      "Dirichlet"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPointNeumann",
      "PointNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroLineNeumann",
      "LineNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> (
      "PoroSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> (
      "PoroVolumeNeumann", "VolumeNeumann"));

  conditions_to_copy.insert(std::pair<std::string, std::string> ("NoPenetration",
      "NoPenetration"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPartInt",
      "PoroPartInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroCoupling",
      "PoroCoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("FSICoupling",
      "FSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("FPSICoupling",
      "FPSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPresInt",
      "PoroPresInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("SurfFlowRate",
      "SurfFlowRate"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("LineFlowRate",
      "LineFlowRate"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(
    const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_fluidporo))
    dserror("Material with ID %d is not admissible for fluid poroelasticity elements",matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element*              oldele,
    const int                  matid,
    const bool                 isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::Fluid> fluid = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Fluid>(newele);
  if (fluid!=Teuchos::null)
  {
    fluid->SetMaterial(matid);
    //Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())->SetInitialPorosity(
              Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->Initporosity());
    fluid->SetDisType(oldele->Shape()); // set distype as well!
    fluid->SetIsAle(true);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support fluid elements here
    eletype.push_back("FLUIDPORO");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support transport elements here
    eletype.push_back("TRANSP");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroScatraCloneStrategy::ConditionsToCopy()
{
  //call base class
  std::map<std::string, std::string> conditions_to_copy = SCATRA::ScatraFluidCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroCoupling",
      "PoroCoupling"));

  return conditions_to_copy;
}



