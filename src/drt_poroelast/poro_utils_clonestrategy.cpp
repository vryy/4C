/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

\maintainer Christoph Ager
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

#include "../drt_so3/so_nurbs27.H"
#include "../drt_so3/so3_scatra.H"
#include "../drt_so3/so3_poro_scatra.H"
#include "../drt_so3/so3_poro_p1_scatra.H"

#include "../drt_w1/wall1_scatra.H"
#include "../drt_w1/wall1_poro_scatra.H"
#include "../drt_w1/wall1_poro_p1_scatra.H"

#include "../drt_s8/shell8_scatra.H"

#include "../drt_scatra_ele/scatra_ele.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroDirichlet", "Dirichlet"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroPointNeumann", "PointNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroLineNeumann", "LineNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroVolumeNeumann", "VolumeNeumann"));

  conditions_to_copy.insert(std::pair<std::string, std::string>("NoPenetration", "NoPenetration"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroPartInt", "PoroPartInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroCoupling", "PoroCoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("FSICoupling", "FSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("FPSICoupling", "FPSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroPresInt", "PoroPresInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("Mortar", "Mortar"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("SurfFlowRate", "SurfFlowRate"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("LineFlowRate", "LineFlowRate"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ImmersedSearchbox", "ImmersedSearchbox"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("XFEMSurfFPIMono", "XFEMSurfFPIMono"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_fluidporo))
    dserror("Material with ID %d is not admissible for fluid poroelasticity elements", matid);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);
  if (fluid != Teuchos::null)
  {
    fluid->SetMaterial(matid);
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())
        ->SetInitialPorosity(
            Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->Initporosity());
    fluid->SetDisType(oldele->Shape());  // set distype as well!
    fluid->SetIsAle(true);
    DRT::ELEMENTS::So_base* so_base = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
    if (so_base)
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
  // clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support fluid elements here
    eletype.push_back("FLUIDPORO");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*
 | return SCATRA::ImplType of element (public)            schmidt 09/17 |
 *----------------------------------------------------------------------*/
INPAR::SCATRA::ImplType POROELAST::UTILS::PoroScatraCloneStrategy::GetImplType(
    DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
)
{
  INPAR::SCATRA::ImplType impltype(INPAR::SCATRA::impltype_undefined);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = ele->ElementType().Name();

  // TET 4 Elements
  // tet 4 solid poro scatra
  if (eletypename == "So_tet4PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(
             ele))
            ->ImplType();
  // tet4 solid porop1 scatra
  else if (eletypename == "So_tet4PoroP1ScatraType")
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(ele))
            ->ImplType();
  // tet 10 solid poro scatra
  else if (eletypename == "So_tet10PoroScatraType")
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>*>(ele))
            ->ImplType();
  // HEX 8 Elements
  // hex8 solid poro scatra
  else if (eletypename == "So_hex8PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
             ele))
            ->ImplType();
  // hex8 solid porop1 scatra
  else if (eletypename == "So_hex8PoroP1ScatraType")
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(ele))
            ->ImplType();
  // hex27 solid poro scatra
  else if (eletypename == "So_hex27PoroScatraType")
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(ele))
            ->ImplType();
  // nurbs 27
  else if (eletypename == "So_nurbs27PoroScatraType")
    impltype = (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
                    DRT::Element::nurbs27>*>(ele))
                   ->ImplType();
  // wall poro scatra elements
  // quad 4
  else if (eletypename == "WallQuad4PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::quad4>*>(ele))->ImplType();
  // quad 9
  else if (eletypename == "WallQuad9PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::quad9>*>(ele))->ImplType();
  // nurbs 4
  else if (eletypename == "WallNurbs4PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::nurbs4>*>(ele))->ImplType();
  // nurbs 9
  else if (eletypename == "WallNurbs9PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::nurbs9>*>(ele))->ImplType();
  // tri 3
  else if (eletypename == "WallTri3PoroScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::tri3>*>(ele))->ImplType();
  // wall poro p1 elements
  // quad 4
  else if (eletypename == "WallQuad4PoroP1ScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>*>(ele))->ImplType();
  // quad 9
  else if (eletypename == "WallQuad9PoroP1ScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>*>(ele))->ImplType();
  // tri 3
  else if (eletypename == "WallTri3PoroP1ScatraType")
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>*>(ele))->ImplType();
  // call base class routine
  else
    impltype = my::GetImplType(ele);

  return impltype;
}


/*----------------------------------------------------------------------*
 |                                                         vuong 08/11  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support transport elements here
    eletype.push_back("TRANSP");
    return true;
  }

  return false;
}


/*----------------------------------------------------------------------*
 | set the element data (protected)                       schmidt 09/17 |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroScatraCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != NULL)
  {
    // set material
    trans->SetMaterial(matid, oldele);
    // set distype as well!
    trans->SetDisType(oldele->Shape());

    // now check whether ImplType is reasonable and if set the ImplType
    INPAR::SCATRA::ImplType impltype =
        POROELAST::UTILS::PoroScatraCloneStrategy::GetImplType(oldele);
    if (impltype == INPAR::SCATRA::impltype_undefined)
      dserror(
          "PoroScatraCloneStrategy copies scatra discretization from structure discretization, but "
          "the "
          "STRUCTURE elements that are defined in the .dat file are either not meant to be copied "
          "to scatra elements "
          "or the ImplType is set 'Undefined' which is not meaningful for the created scatra "
          "discretization! "
          "Use SOLIDSCATRA, WALLSCATRA, SHELLSCATRA, SOLIDPOROSCATRA, SOLIDPOROP1SCATRA, "
          "WALLPOROSCATRA or "
          "WALLPOROP1SCATRA Elements with meaningful ImplType instead!");

    else
      trans->SetImplType(impltype);
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
std::map<std::string, std::string> POROELAST::UTILS::PoroScatraCloneStrategy::ConditionsToCopy()
{
  // call base class
  std::map<std::string, std::string> conditions_to_copy =
      SSI::ScatraStructureCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroCoupling", "PoroCoupling"));

  conditions_to_copy.insert(std::pair<std::string, std::string>("Initfield", "Initfield"));

  // artery to scatra coupling
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ArtScatraCouplCon", "ArtScatraCouplCon"));

  // oxygen partial pressure calculation
  conditions_to_copy.insert(std::pair<std::string, std::string>(
      "PoroMultiphaseScatraOxyPartPressCalcCond", "PoroMultiphaseScatraOxyPartPressCalcCond"));

  // Robin boundary condition
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportRobin", "TransportRobin"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*
 |                                                         rauch 03/15  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastImmersedCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
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
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoroImmersed> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed>(newele);
  if (fluid != Teuchos::null)
  {
    fluid->SetMaterial(matid);
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())
        ->SetInitialPorosity(
            Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->Initporosity());
    fluid->SetDisType(oldele->Shape());  // set distype as well!
    fluid->SetIsAle(true);
    DRT::ELEMENTS::So_base* so_base = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
    if (so_base)
    {
      fluid->SetKinematicType(so_base->KinematicType());
      if (so_base->KinematicType() == INPAR::STR::kinem_vague) dserror("undefined kinematic type");
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
