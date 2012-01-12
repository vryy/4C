/*!------------------------------------------------------------------------------------------------*
\file topopt_utils.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_UTILS_CPP_
#define TOPOPT_UTILS_CPP_

#include <iostream>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_scatra/scatra_element.H"
#include "topopt_utils.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<string,string> TOPOPT::TopoptFluidCloneStrategy::ConditionsToCopy()
{
  std::map<string,string> conditions_to_copy;
// TODO adopt these conditions if optimization field is no more equivalent from fluid (winklmaier)
  conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<string,string>("TransportPointNeumann","PointNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportLineNeumann","LineNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportSurfaceNeumann","SurfaceNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportVolumeNeumann","VolumeNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportNeumannInflow","TransportNeumannInflow"));
  // when the fluid problem is periodic we also expect the mass transport to be so:
  conditions_to_copy.insert(pair<string,string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));

  conditions_to_copy.insert(pair<string,string>("LineNeumann","FluidLineNeumann"));
  conditions_to_copy.insert(pair<string,string>("SurfaceNeumann","FluidSurfaceNeumann"));
  conditions_to_copy.insert(pair<string,string>("VolumeNeumann","FluidVolumeNeumann"));
  conditions_to_copy.insert(pair<string,string>("KrylovSpaceProjection","KrylovSpaceProjection"));
  conditions_to_copy.insert(pair<string,string>("Initfield","Initfield"));

  // for moving boundary problems
  conditions_to_copy.insert(pair<string,string>("FSICoupling","FSICoupling"));

  // meshtying
  conditions_to_copy.insert(pair<string,string>("Contact","Contact"));

  return conditions_to_copy;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TOPOPT::TopoptFluidCloneStrategy::CheckMaterialType(const int matid)
{
// We take the material with the ID specified by the user
// Here we check first, whether this material is of admissible type
INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();

if (mtype != INPAR::MAT::m_opti_dens)
{
  dserror("Material with ID %d is not admissible for optimization density elements",matid);
}
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TOPOPT::TopoptFluidCloneStrategy::SetElementData(
    RCP<DRT::Element> newele,
    DRT::Element* oldele,
    const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
#if defined(D_FLUID3)
      //TODO fix this hack -> own element class for topopt? or fix this at higher level... (winklmaier)
      //      if fixed remove the scatra include
      DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
      if (trans!=NULL)
      {
        trans->SetMaterial(matid);
        trans->SetDisType(oldele->Shape()); // set distype as well!
      }
      else
#endif
    {
      dserror("unsupported element type '%s'", typeid(*newele).name());
    }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool TOPOPT::TopoptFluidCloneStrategy::DetermineEleType(
    DRT::Element* actele,
    const bool ismyele,
    vector<string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true; // yes, we copy EVERY element (no submeshes)
}



void TOPOPT::printTopOptLogo()
{
  std::cout << "      _________            " << std::endl;
  std::cout << "     /         \\          " << std::endl;
  std::cout << "    /   _____   \\         " << " Das ist das       " << std::endl;
  std::cout << "   |   /     \\@ @\\        " << " Gebiets-          " << std::endl;
  std::cout << "   |   |      \\__/        " << " Optimierungsmodul " << std::endl;
  std::cout << "   |   \\         \\        " << " in BACI           " << std::endl;
  std::cout << "    \\   \\         \\        " << "                   " << std::endl;
  std::cout << "     \\   \\                " << " Die Schlange      " << std::endl;
  std::cout << "      \\   \\_________      " << " wird sich bald    " << std::endl;
  std::cout << "       \\             \\    " << " teilen und wieder " << std::endl;
  std::cout << "        \\__________   \\   " << " zusammenwachsen   " << std::endl;
  std::cout << "                   \\  |   " << " können!           " << std::endl;
  std::cout << "        _          |  |    " << std::endl;
  std::cout << "       | \\        /  /    " << std::endl;
  std::cout << "        \\ \\______/  /     " << std::endl;
  std::cout << "         \\_________/      " << std::endl;
  std::cout << "                           " << std::endl;
}

#endif /* TOPOPT_UTILS_CPP_ */
#endif  // #ifdef CCADISCRET
