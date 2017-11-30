/*!
\file scatra_ele_calc_utils.cpp

\brief Utility methods for scatra

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/


#include "scatra_ele_calc_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_condition_utils.H"

namespace SCATRA
{

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool IsBinaryElectrolyte(const std::vector<double>& valence)
{
  int numions(0);
  for (size_t k=0; k < valence.size(); k++)
  {
    if (abs(valence[k]) > EPS10)
      numions++;
  }
  return (numions == 2);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> GetIndicesBinaryElectrolyte(const std::vector<double>& valence)
{
  // indices of the two charged species to be determined
  std::vector<int> indices;
  for (size_t k=0; k < valence.size(); k++)
  {
    // is there some charge?
    if (abs(valence[k]) > EPS10)
      indices.push_back(k);
  }
  if (indices.size() != 2) dserror("Found no binary electrolyte!");

  return indices;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double CalResDiffCoeff(
    const std::vector<double>& valence,
    const std::vector<double>& diffus,
    const std::vector<int>& indices
    )
{
  if (indices.size() != 2) dserror("Non-matching number of indices!");
  const int first = indices[0];
  const int second = indices[1];
  if ((valence[first]*valence[second])>EPS10)
    dserror("Binary electrolyte has no opposite charges.");
  const double n = ((diffus[first]*valence[first])-(diffus[second]*valence[second]));
  if (abs(n) < EPS12)
    dserror("denominator in resulting diffusion coefficient is nearly zero");

  return diffus[first]*diffus[second]*(valence[first]-valence[second])/n;
}


/*-------------------------------------------------------------------------------*
 |find elements of inflow section                                rasthofer 01/12 |
 |for turbulent low Mach number flows with turbulent inflow condition            |
 *-------------------------------------------------------------------------------*/
bool InflowElement(const DRT::Element* ele)
{
  bool inflow_ele = false;

  std::vector<DRT::Condition*> myinflowcond;

  // check whether all nodes have a unique inflow condition
  DRT::UTILS::FindElementConditions(ele, "TurbulentInflowSection", myinflowcond);
  if (myinflowcond.size()>1)
    dserror("More than one inflow condition on one node!");

  if (myinflowcond.size()==1)
   inflow_ele = true;

  return inflow_ele;
}


/*---------------------------------------------------------------------------------------------------------------------*
 | convert implementation type of scalar transport elements into corresponding string for output purposes   fang 02/15 |
 *---------------------------------------------------------------------------------------------------------------------*/
const std::string ImplTypeToString(const INPAR::SCATRA::ImplType impltype)
{
  // determine implementation type
  std::string impltypestring;
  switch(impltype)
  {
  case INPAR::SCATRA::impltype_std:
  {
    impltypestring = "Standard scalar transport";
    break;
  }
  case INPAR::SCATRA::impltype_thermo_elch_electrode:
  {
    impltypestring = "Heat transport within electrodes";
    break;
  }
  case INPAR::SCATRA::impltype_thermo_elch_diffcond:
  {
    impltypestring = "Heat transport within concentrated electrolytes";
    break;
  }
  case INPAR::SCATRA::impltype_advreac:
  {
    impltypestring = "Advanced reactions";
    break;
  }
  case INPAR::SCATRA::impltype_refconcreac:
  {
    impltypestring = "Reference concentrations AND reactions";
    break;
  }
  case INPAR::SCATRA::impltype_chemo:
  {
     impltypestring = "Chemotaxis";
     break;
  }
  case INPAR::SCATRA::impltype_chemoreac:
  {
     impltypestring = "Advanced reactions AND chemotaxis";
     break;
  }
  case INPAR::SCATRA::impltype_aniso:
  {
    impltypestring = "Anisotropic scalar transport";
    break;
  }
  case INPAR::SCATRA::impltype_cardiac_monodomain:
  {
    impltypestring = "Cardiac monodomain";
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond:
  {
    impltypestring = "Electrochemistry for diffusion-conduction formulation";
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond_multiscale:
  {
    impltypestring = "Electrochemistry for diffusion-conduction formulation within a multi-scale framework";
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond_thermo:
  {
    impltypestring = "Electrochemistry for diffusion-conduction formulation with thermal effects";
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode:
  {
    impltypestring = "Electrochemistry for electrodes";
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode_growth:
  {
    impltypestring = "Electrochemistry for electrodes exhibiting lithium plating and stripping";
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode_thermo:
  {
    impltypestring = "Electrochemistry for electrodes with thermal effects";
    break;
  }
  case INPAR::SCATRA::impltype_elch_NP:
  {
    impltypestring = "Electrochemistry for Nernst-Planck formulation";
    break;
  }
  case INPAR::SCATRA::impltype_loma:
  {
    impltypestring = "Low Mach number flow";
    break;
  }
  case INPAR::SCATRA::impltype_levelset:
  {
    impltypestring = "Levelset without reinitialization";
    break;
  }
  case INPAR::SCATRA::impltype_lsreinit:
  {
    impltypestring = "Levelset with reinitialization";
    break;
  }
  case INPAR::SCATRA::impltype_poro:
  {
    impltypestring = "Scalar transport in porous media";
    break;
  }
  case INPAR::SCATRA::impltype_pororeac:
  {
    impltypestring = "Reactive scalar transport in porous media";
    break;
  }
  case INPAR::SCATRA::impltype_bondreac:
  {
    impltypestring = "Advanced reactions with bond dynamics";
    break;
  }
  case INPAR::SCATRA::impltype_variational_diffusion:
  {
    impltypestring = "Chemical diffusion under a variational setting";
    break;
  }
  case INPAR::SCATRA::impltype_undefined:
  {
    impltypestring = "Undefined";
    break;
  }
  default:
  {
    dserror("Invalid implementation type!");
    break;
  }
  }

  return impltypestring;
}
} // namespace SCATRA
