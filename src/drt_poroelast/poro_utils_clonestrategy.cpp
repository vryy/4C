/*----------------------------------------------------------------------*/
/*!
 \file poro_utils_clonestrategy.cpp

 \brief strategy to clone porofluid from porous solid

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_utils_clonestrategy.H"

#include "poroelast_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid_ele/fluid_ele_poro_immersed.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include"../drt_w1/wall1_poro.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
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
  conditions_to_copy.insert(std::pair<std::string, std::string> ("Mortar",
        "Mortar"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("SurfFlowRate",
      "SurfFlowRate"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("LineFlowRate",
      "LineFlowRate"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("ImmersedSearchbox",
      "ImmersedSearchbox"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("XFEMSurfFPIMono",
      "XFEMSurfFPIMono"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(
    const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_fluidporo))
    dserror("Material with ID %d is not admissible for fluid poroelasticity elements",matid);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element*              oldele,
    const int                  matid,
    const bool                 isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);
  if (fluid!=Teuchos::null)
  {
    fluid->SetMaterial(matid);
    //Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())->SetInitialPorosity(
              Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->Initporosity());
    fluid->SetDisType(oldele->Shape()); // set distype as well!
    fluid->SetIsAle(true);
    DRT::ELEMENTS::So_base*  so_base  = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
    if(so_base)
      fluid->SetKinematicType(so_base->KinematicType());
    else
      dserror(" dynamic cast from DRT::Element* to DRT::ELEMENTS::So_base* failed ");
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroScatraCloneStrategy::ConditionsToCopy()
{
  //call base class
  std::map<std::string, std::string> conditions_to_copy = SCATRA::ScatraFluidCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroCoupling",
      "PoroCoupling"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*
 |                                                         rauch 03/15  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastImmersedCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support transport elements here
    eletype.push_back("FLUIDPOROIMMERSED");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*
 |                                                         rauch 03/15  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastImmersedCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element*              oldele,
    const int                  matid,
    const bool                 isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoroImmersed> fluid = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed>(newele);
  if (fluid!=Teuchos::null)
  {
    fluid->SetMaterial(matid);
    //Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())->SetInitialPorosity(
              Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->Initporosity());
    fluid->SetDisType(oldele->Shape()); // set distype as well!
    fluid->SetIsAle(true);
    DRT::ELEMENTS::So_base*  so_base  = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
    if(so_base)
    {
      fluid->SetKinematicType(so_base->KinematicType());
      if(so_base->KinematicType()==INPAR::STR::kinem_vague)
        dserror("undefined kinematic type");
    }
    else
      dserror(" dynamic cast from DRT::Element* to DRT::ELEMENTS::So_base* failed ");
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}
