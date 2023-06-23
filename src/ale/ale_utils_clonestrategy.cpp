/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone ALE discretization form other discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale_utils_clonestrategy.H"
#include "lib_globalproblem.H"
#include "mat_par_material.H"
#include "mat_par_bundle.H"

// we need to know all element types for the ale mesh creation
#include "fluid_ele.H"
#include "ale_ale2.H"
#include "ale_ale2_nurbs.H"
#include "ale_ale3.H"
#include "ale_ale3_nurbs.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::map<std::string, std::string> ALE::UTILS::AleCloneStrategy::ConditionsToCopy() const
{
  return {{"ALEDirichlet", "Dirichlet"}, {"FSICoupling", "FSICoupling"},
      {"FPSICoupling", "FPSICoupling"}, {"FREESURFCoupling", "FREESURFCoupling"},
      {"ALEUPDATECoupling", "ALEUPDATECoupling"}, {"StructAleCoupling", "StructAleCoupling"},
      {"LinePeriodic", "LinePeriodic"}, {"SurfacePeriodic", "SurfacePeriodic"},
      {"ElchBoundaryKinetics", "ElchBoundaryKinetics"},
      {"XFEMSurfFluidFluid", "XFEMSurfFluidFluid"}, {"FluidFluidCoupling", "FluidFluidCoupling"},
      {"AleWear", "AleWear"}, {"AleLocsys", "Locsys"}, {"Mortar", "Mortar"},
      {"UncertainSurface", "UncertainSurface"}};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::AleCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if (mtype != INPAR::MAT::m_stvenant && mtype != INPAR::MAT::m_elasthyper)
    dserror("Material with ID %d is not admissible for ALE elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::AleCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool nurbsdis)
{
  if (nurbsdis == false)
  {
    DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
    if (ale2 != NULL)
    {
      ale2->SetMaterial(matid);
    }
    else
    {
      DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(newele.get());
      if (ale3 != NULL)
      {
        ale3->SetMaterial(matid);
      }
      else
      {
        dserror("unsupported ale element type '%s'", typeid(*newele).name());
      }
    }
  }
  else
  {
    DRT::ELEMENTS::NURBS::Ale2Nurbs* ale2 =
        dynamic_cast<DRT::ELEMENTS::NURBS::Ale2Nurbs*>(newele.get());
    if (ale2 != NULL)
    {
      ale2->SetMaterial(matid);
    }
    else
    {
      DRT::ELEMENTS::NURBS::Ale3Nurbs* ale3 =
          dynamic_cast<DRT::ELEMENTS::NURBS::Ale3Nurbs*>(newele.get());

      if (ale3 != NULL)
      {
        ale3->SetMaterial(matid);
      }
      else
      {
        dserror("unsupported ale element type '%s'", typeid(*newele).name());
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool ALE::UTILS::AleCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  bool cloneit = true;

  // Fluid meshes may be split into Eulerian and ALE regions.
  // Check, whether actele is a fluid element in order to account for
  // the possible split in Eulerian an ALE regions
  DRT::ELEMENTS::Fluid* f3 = dynamic_cast<DRT::ELEMENTS::Fluid*>(actele);
  if (f3 != NULL)
  {
    cloneit = f3->IsAle();  // if not ALE, element will not be cloned
                            // --> theoretically, support of Eulerian sub meshes
  }

  // Clone it now.
  if (cloneit and ismyele)
  {
    const int nsd = CORE::DRT::UTILS::getDimension(actele->Shape());
    if (nsd == 3)
      eletype.push_back("ALE3");
    else if (nsd == 2)
      eletype.push_back("ALE2");
    else
      dserror("%i D Dimension not supported", nsd);
  }

  return cloneit;
}
