/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone ALE discretization form other discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_utils_clonestrategy.hpp"

#include "4C_ale_ale2.hpp"
#include "4C_ale_ale2_nurbs.hpp"
#include "4C_ale_ale3.hpp"
#include "4C_ale_ale3_nurbs.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::map<std::string, std::string> ALE::UTILS::AleCloneStrategy::conditions_to_copy() const
{
  return {{"ALEDirichlet", "Dirichlet"}, {"FSICoupling", "FSICoupling"},
      {"fpsi_coupling", "fpsi_coupling"}, {"FREESURFCoupling", "FREESURFCoupling"},
      {"ALEUPDATECoupling", "ALEUPDATECoupling"}, {"StructAleCoupling", "StructAleCoupling"},
      {"LinePeriodic", "LinePeriodic"}, {"SurfacePeriodic", "SurfacePeriodic"},
      {"ElchBoundaryKinetics", "ElchBoundaryKinetics"},
      {"XFEMSurfFluidFluid", "XFEMSurfFluidFluid"}, {"FluidFluidCoupling", "FluidFluidCoupling"},
      {"AleWear", "AleWear"}, {"AleLocsys", "Locsys"}, {"Mortar", "Mortar"},
      {"UncertainSurface", "UncertainSurface"}};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::AleCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  CORE::Materials::MaterialType mtype =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if (mtype != CORE::Materials::m_stvenant && mtype != CORE::Materials::m_elasthyper)
    FOUR_C_THROW("Material with ID %d is not admissible for ALE elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::AleCloneStrategy::set_element_data(Teuchos::RCP<CORE::Elements::Element> newele,
    CORE::Elements::Element* oldele, const int matid, const bool nurbsdis)
{
  if (nurbsdis == false)
  {
    DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
    if (ale2 != nullptr)
    {
      ale2->SetMaterial(0, MAT::Factory(matid));
    }
    else
    {
      DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(newele.get());
      if (ale3 != nullptr)
      {
        ale3->SetMaterial(0, MAT::Factory(matid));
      }
      else
      {
        FOUR_C_THROW("unsupported ale element type '%s'", typeid(*newele).name());
      }
    }
  }
  else
  {
    DRT::ELEMENTS::NURBS::Ale2Nurbs* ale2 =
        dynamic_cast<DRT::ELEMENTS::NURBS::Ale2Nurbs*>(newele.get());
    if (ale2 != nullptr)
    {
      ale2->SetMaterial(0, MAT::Factory(matid));
    }
    else
    {
      DRT::ELEMENTS::NURBS::Ale3Nurbs* ale3 =
          dynamic_cast<DRT::ELEMENTS::NURBS::Ale3Nurbs*>(newele.get());

      if (ale3 != nullptr)
      {
        ale3->SetMaterial(0, MAT::Factory(matid));
      }
      else
      {
        FOUR_C_THROW("unsupported ale element type '%s'", typeid(*newele).name());
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool ALE::UTILS::AleCloneStrategy::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  bool cloneit = true;

  // Fluid meshes may be split into Eulerian and ALE regions.
  // Check, whether actele is a fluid element in order to account for
  // the possible split in Eulerian an ALE regions
  DRT::ELEMENTS::Fluid* f3 = dynamic_cast<DRT::ELEMENTS::Fluid*>(actele);
  if (f3 != nullptr)
  {
    cloneit = f3->IsAle();  // if not ALE, element will not be cloned
                            // --> theoretically, support of Eulerian sub meshes
  }

  // Clone it now.
  if (cloneit and ismyele)
  {
    const int nsd = CORE::FE::getDimension(actele->Shape());
    if (nsd == 3)
      eletype.push_back("ALE3");
    else if (nsd == 2)
      eletype.push_back("ALE2");
    else
      FOUR_C_THROW("%i D Dimension not supported", nsd);
  }

  return cloneit;
}

FOUR_C_NAMESPACE_CLOSE
