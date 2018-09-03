/*----------------------------------------------------------------------------*/
/*!
\file ale_utils_clonestrategy.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale_utils_clonestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"

// we need to know all element types for the ale mesh creation
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_ale2/ale2.H"
#include "../drt_ale2/ale2_nurbs.H"
#include "../drt_ale3/ale3.H"
#include "../drt_ale3/ale3_nurbs.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::map<std::string, std::string> ALE::UTILS::AleCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string>("ALEDirichlet", "Dirichlet"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("FSICoupling", "FSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("FPSICoupling", "FPSICoupling"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("FREESURFCoupling", "FREESURFCoupling"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ALEUPDATECoupling", "ALEUPDATECoupling"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("StructAleCoupling", "StructAleCoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("LinePeriodic", "LinePeriodic"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("SurfacePeriodic", "SurfacePeriodic"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ElchBoundaryKinetics", "ElchBoundaryKinetics"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("XFEMSurfFluidFluid", "XFEMSurfFluidFluid"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("FluidFluidCoupling", "FluidFluidCoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("AleWear", "AleWear"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("AleLocsys", "Locsys"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("Mortar", "Mortar"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("UncertainSurface", "UncertainSurface"));

  return conditions_to_copy;
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
#ifdef D_ALE
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
#endif

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
    const int nsd = DRT::UTILS::getDimension(actele->Shape());
    if (nsd == 3)
      eletype.push_back("ALE3");
    else if (nsd == 2)
      eletype.push_back("ALE2");
    else
      dserror("%i D Dimension not supported", nsd);
  }

  return cloneit;
}
