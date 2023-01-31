/*----------------------------------------------------------------------*/
/*! \file
\brief mesh clone strategy for scatra-structure-interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/

#include "ssi_clonestrategy.H"

#include "inpar_ssi.H"

#include "lib_globalproblem.H"
#include "lib_discret.H"

#include "mat_par_material.H"
#include "mat_par_bundle.H"

#include "scatra_ele.H"

#include "adapter_structure_scatra_ele.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> SSI::ScatraStructureCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string>("TransportDirichlet", "Dirichlet"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportPointNeumann", "PointNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportLineNeumann", "LineNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportVolumeNeumann", "VolumeNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TransportNeumannInflow", "TransportNeumannInflow"));
  // for moving boundary problems
  conditions_to_copy.insert(std::pair<std::string, std::string>("FSICoupling", "FSICoupling"));
  // for coupled scalar transport fields
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ScaTraCoupling", "ScaTraCoupling"));
  // boundary flux evaluation condition for scalar transport
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ScaTraFluxCalc", "ScaTraFluxCalc"));
  // initial field conditions
  conditions_to_copy.insert(std::pair<std::string, std::string>("Initfield", "Initfield"));
  // copy s2i mesh tying condition from its ssi equivalent
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("SSIInterfaceMeshtying", "S2IMeshtying"));
  // scatra-scatra interface kinetics conditions
  conditions_to_copy.insert(std::pair<std::string, std::string>("S2IKinetics", "S2IKinetics"));
  // copy partitioning of the scatra field
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ScatraPartitioning", "ScatraPartitioning"));
  // electrode state of charge conditions
  conditions_to_copy.insert(std::pair<std::string, std::string>("ElectrodeSOC", "ElectrodeSOC"));
  // cell voltage conditions
  conditions_to_copy.insert(std::pair<std::string, std::string>("CellVoltage", "CellVoltage"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("CellVoltagePoint", "CellVoltagePoint"));
  // calculate total and mean scalar
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("TotalAndMeanScalar", "TotalAndMeanScalar"));
  // cell cycling conditions
  conditions_to_copy.insert(std::pair<std::string, std::string>("CCCVCycling", "CCCVCycling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("CCCVHalfCycle", "CCCVHalfCycle"));
  // convective heat transfer conditions (Newton's law of heat transfer)
  conditions_to_copy.insert(std::pair<std::string, std::string>(
      "TransportThermoConvections", "TransportThermoConvections"));
  conditions_to_copy.insert(std::pair<std::string, std::string>(
      "SSIMeshtying3DomainIntersection", "Meshtying3DomainIntersection"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("SSISurfaceManifold", "SSISurfaceManifold"));
  conditions_to_copy.insert(std::pair<std::string, std::string>(
      "SSISurfaceManifoldKinetics", "SSISurfaceManifoldKinetics"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> SSI::ScatraStructureCloneStrategyManifold::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;
  return conditions_to_copy;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
INPAR::SCATRA::ImplType SSI::ScatraStructureCloneStrategy::GetImplType(DRT::Element* ele)
{
  return ADAPTER::GetScaTraImplType(ele);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_elchmat) &&
      (mtype != INPAR::MAT::m_electrode) && (mtype != INPAR::MAT::m_matlist) &&
      (mtype != INPAR::MAT::m_matlist_reactions) && (mtype != INPAR::MAT::m_myocard) &&
      (mtype != INPAR::MAT::m_thermostvenant))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set material
    trans->SetMaterial(matid, oldele);
    // set distype as well
    trans->SetDisType(oldele->Shape());

    // now check whether ImplType is reasonable and if set the ImplType
    INPAR::SCATRA::ImplType impltype = SSI::ScatraStructureCloneStrategy::GetImplType(oldele);

    if (impltype == INPAR::SCATRA::impltype_undefined)
    {
      dserror(
          "ScatraStructureCloneStrategy copies scatra discretization from structure "
          "discretization, but the STRUCTURE elements that are defined in the .dat file are "
          "either not meant to be copied to scatra elements or the ImplType is set 'Undefined' "
          "which is not meaningful for the created scatra discretization! Use SOLIDSCATRA, "
          "WALLSCATRA, SHELLSCATRA or TRUSS3SCATRA elements with meaningful ImplType instead!");
    }
    else
      trans->SetImplType(impltype);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureCloneStrategyManifold::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // determine impl type from manifold condition by identifying the condition for this element
  auto struct_dis = DRT::Problem::Instance()->GetDis("structure");

  std::vector<DRT::Condition*> conditions;
  struct_dis->GetCondition("SSISurfaceManifold", conditions);

  auto impltype = INPAR::SCATRA::impltype_undefined;
  for (auto* condition : conditions)
  {
    auto cond_eles = condition->Geometry();
    if (cond_eles.find(oldele->Id()) != cond_eles.end())
    {
      impltype = static_cast<INPAR::SCATRA::ImplType>(condition->GetInt("ImplType"));
      continue;
    }
  }

  if (impltype != INPAR::SCATRA::impltype_elch_electrode and
      impltype != INPAR::SCATRA::impltype_elch_diffcond and impltype != INPAR::SCATRA::impltype_std)
    dserror("Scatra Impltype not supported for SSI with transport on manifolds");

  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans == nullptr or oldele->ElementType().Name() != "StructuralSurfaceType")
    dserror("element type not supported");

  // set material
  trans->SetMaterial(matid, oldele);
  // set distype
  trans->SetDisType(oldele->Shape());
  // set impltype
  trans->SetImplType(impltype);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::ScatraStructureCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.emplace_back("TRANSP");

  return true;  // yes, we copy EVERY element (no submeshes)
}
