/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for cloning the porofluid discretization

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_poromultiphase_utils_clonestrategy.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_porofluidmultiphase_ele.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | define conditions to copy to the cloned discretization    vuong 08/16 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string>
POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::conditions_to_copy() const
{
  return {{"PoroDirichlet", "Dirichlet"}, {"PoroPointNeumann", "PointNeumann"},
      {"PoroLineNeumann", "LineNeumann"}, {"PoroSurfaceNeumann", "SurfaceNeumann"},
      {"PoroVolumeNeumann", "VolumeNeumann"}, {"Initfield", "Initfield"},
      {"ArtPorofluidCouplConNodebased", "ArtPorofluidCouplConNodebased"},
      {"ArtPorofluidCouplConNodeToPoint", "ArtPorofluidCouplConNodeToPoint"}};
}


/*----------------------------------------------------------------------*
 | check for correct material                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if ((mtype != Core::Materials::m_fluidporo_multiphase) and
      (mtype != Core::Materials::m_fluidporo_multiphase_reactions))
    FOUR_C_THROW("Material with ID %d is not admissible for porofluid multiphase elements", matid);
}


/*----------------------------------------------------------------------*
 | set element-specific data (material etc.)                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::set_element_data(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele, const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the PoroFluidMultiPhase element!
  Discret::ELEMENTS::PoroFluidMultiPhase* porofluidele =
      dynamic_cast<Discret::ELEMENTS::PoroFluidMultiPhase*>(newele.get());
  if (porofluidele != nullptr)
  {
    porofluidele->SetMaterial(0, Mat::Factory(matid));
    porofluidele->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*
 | determine whether element is copied or not               vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support POROFLUIDMULTIPHASE elements here
  eletype.push_back("POROFLUIDMULTIPHASE");

  // all elements are copied
  return true;
}

FOUR_C_NAMESPACE_CLOSE
