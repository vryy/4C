/*!------------------------------------------------------------------------------------------------*
\file topopt_utils.cpp

\brief collection of functions in namespace topopt

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_utils.H"
#include "topopt_optimizer_ele.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
map<string,string> TOPOPT::TopoptFluidCloneStrategy::ConditionsToCopy()
{
  map<string,string> conditions_to_copy;

  // when the fluid problem is periodic we also expect the optimization to be so:
  conditions_to_copy.insert(std::pair<string,string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(std::pair<string,string>("SurfacePeriodic","SurfacePeriodic"));

  // parts of the objective might use inflow and outflow integrals. Therefore, we
  // save potential inflow and outflow regions as conditions for the optimization

  // neumann conditions are supposed to be outflow regions (TODO: handle neumann inflow)
  conditions_to_copy.insert(std::pair<string,string>("LineNeumann","LineOutflow"));
  conditions_to_copy.insert(std::pair<string,string>("SurfaceNeumann","SurfaceOutflow"));

  // dirichlet conditions are supposed to be inflow regions
  conditions_to_copy.insert(std::pair<string,string>("LineDirichlet","LineInflow"));
  conditions_to_copy.insert(std::pair<string,string>("SurfaceNeumann","SurfaceInflow"));

  // initial field also in optimization possibly
  conditions_to_copy.insert(std::pair<string,string>("Initfield","Initfield"));

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
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element* oldele,
    const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  DRT::ELEMENTS::TopOpt* trans = dynamic_cast<DRT::ELEMENTS::TopOpt*>(newele.get());
  if (trans!=NULL)
  {
    trans->SetMaterial(matid);
    trans->SetDisType(oldele->Shape()); // set distype as well!
  }
  else
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
    std::vector<string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TOPOPT");

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
  std::cout << "                   \\  |   " << " koennen!           " << std::endl;
  std::cout << "        _          |  |    " << std::endl;
  std::cout << "       | \\        /  /    " << std::endl;
  std::cout << "        \\ \\______/  /     " << std::endl;
  std::cout << "         \\_________/      " << std::endl;
  std::cout << "                           " << std::endl;
}
