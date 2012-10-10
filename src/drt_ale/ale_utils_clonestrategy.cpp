/*----------------------------------------------------------------------*/
/*!
\file ale_utils_clonestrategy.cpp

\brief mesh clone strategy for ale problems

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/


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



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<string,string> ALE::UTILS::AleFluidCloneStrategy::ConditionsToCopy()
{
  std::map<string,string> conditions_to_copy;

  conditions_to_copy.insert(pair<string,string>("ALEDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<string,string>("FSICoupling","FSICoupling"));
  conditions_to_copy.insert(pair<string,string>("FREESURFCoupling","FREESURFCoupling"));
  conditions_to_copy.insert(pair<string,string>("StructAleCoupling","StructAleCoupling"));
  conditions_to_copy.insert(pair<string,string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));
  conditions_to_copy.insert(pair<string,string>("ElectrodeKinetics","ElectrodeKinetics"));
  conditions_to_copy.insert(pair<string,string>("XFEMCoupling","XFEMCoupling"));
  conditions_to_copy.insert(pair<string,string>("FluidFluidCoupling","FluidFluidCoupling"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::AleFluidCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if (mtype != INPAR::MAT::m_stvenant)
    dserror("Material with ID %d is not admissible for ALE elements",matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::AleFluidCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element* oldele,
    const int matid,
    const bool nurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property.

#ifdef D_ALE
    if(nurbsdis==false)
    {
      DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matid);
      }
      else
      {
        DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(newele.get());
        if (ale3!=NULL)
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
      DRT::ELEMENTS::NURBS::Ale2Nurbs* ale2 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale2Nurbs*>(newele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matid);
      }
      else
      {
        DRT::ELEMENTS::NURBS::Ale3Nurbs* ale3 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale3Nurbs*>(newele.get());

        if(ale3!=NULL)
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool ALE::UTILS::AleFluidCloneStrategy::DetermineEleType(
    DRT::Element* actele,
    const bool ismyele,
    vector<string>& eletype)
{
  bool isale = false;
  bool found = false;

#ifdef D_FLUID3
    DRT::ELEMENTS::Fluid* f3 = dynamic_cast<DRT::ELEMENTS::Fluid*>(actele);
    if (not found and f3!=NULL)
    {
      const int  nsd = DRT::UTILS::getDimension(f3->Shape());
      found = true;
      isale = f3->IsAle();

      if (isale and ismyele)
      {
        if (nsd == 3)
          eletype.push_back("ALE3");
        else if (nsd == 2)
          eletype.push_back("ALE2");
        else
          dserror("%i D Dimension not supported", nsd);
      }
    }
#endif

    if (not found)
      dserror("unsupported fluid element type '%s'", typeid(*actele).name());

  return isale;
}

