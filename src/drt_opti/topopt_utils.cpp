/*---------------------------------------------------------------------*/
/*!

\brief collection of functions in namespace topopt

\level 3

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/


#include "topopt_utils.H"
#include "topopt_optimizer_ele.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> TOPOPT::TopoptFluidCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  // when the fluid problem is periodic we also expect the optimization to be so:
  conditions_to_copy.insert(std::pair<std::string, std::string>("LinePeriodic", "LinePeriodic"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("SurfacePeriodic", "SurfacePeriodic"));

  // parts of the objective might use inflow and outflow integrals. Therefore, we
  // save potential inflow and outflow regions as conditions for the optimization

  // neumann conditions are supposed to be outflow regions (TODO: handle neumann inflow)
  conditions_to_copy.insert(std::pair<std::string, std::string>("LineNeumann", "LineOutflow"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("SurfaceNeumann", "SurfaceOutflow"));

  // dirichlet conditions are supposed to be inflow regions
  conditions_to_copy.insert(std::pair<std::string, std::string>("LineDirichlet", "LineInflow"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("SurfaceNeumann", "SurfaceInflow"));

  // initial field also in optimization possibly
  conditions_to_copy.insert(std::pair<std::string, std::string>("Initfield", "Initfield"));

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
    dserror("Material with ID %d is not admissible for optimization density elements", matid);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TOPOPT::TopoptFluidCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  DRT::ELEMENTS::TopOpt* trans = dynamic_cast<DRT::ELEMENTS::TopOpt*>(newele.get());
  if (trans != NULL)
  {
    trans->SetMaterial(matid);
    trans->SetDisType(oldele->Shape());  // set distype as well!
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
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TOPOPT");

  return true;  // yes, we copy EVERY element (no submeshes)
}



void TOPOPT::printTopOptLogo()
{
  std::cout << "      _________            " << std::endl;
  std::cout << "     /         \\          " << std::endl;
  std::cout << "    /   _____   \\         "
            << " Das ist das       " << std::endl;
  std::cout << "   |   /     \\@ @\\        "
            << " Gebiets-          " << std::endl;
  std::cout << "   |   |      \\__/        "
            << " Optimierungsmodul " << std::endl;
  std::cout << "   |   \\         \\        "
            << " in BACI           " << std::endl;
  std::cout << "    \\   \\         \\        "
            << "                   " << std::endl;
  std::cout << "     \\   \\                "
            << " Die Schlange      " << std::endl;
  std::cout << "      \\   \\_________      "
            << " wird sich bald    " << std::endl;
  std::cout << "       \\             \\    "
            << " teilen und wieder " << std::endl;
  std::cout << "        \\__________   \\   "
            << " zusammenwachsen   " << std::endl;
  std::cout << "                   \\  |   "
            << " koennen!           " << std::endl;
  std::cout << "        _          |  |    " << std::endl;
  std::cout << "       | \\        /  /    " << std::endl;
  std::cout << "        \\ \\______/  /     " << std::endl;
  std::cout << "         \\_________/      " << std::endl;
  std::cout << "                           " << std::endl;
}



/// expands a filename
std::string TOPOPT::modifyFilename(
    const std::string& filename, const std::string& expansion, bool handleRestart, bool firstcall)
{
  size_t pos = filename.rfind('/');
  std::string filenameout = filename.substr(0, pos + 1) + expansion + filename.substr(pos + 1);
  if (handleRestart)
  {
    pos = filenameout.rfind('-');
    if (pos != std::string::npos)
    {
      int number = 0;
      filenameout = filenameout.substr(0, pos);
      for (;;)
      {
        number += 1;
        std::stringstream helpname;
        helpname << filenameout << "-" << number << ".control";
        std::ifstream file(helpname.str().c_str());
        if (not file)
        {
          std::stringstream name;
          if (firstcall)
            name << filenameout << "-" << number - 1 << ".control";  // at restart -1 required here
          else
            name << filenameout << "-" << number - 2 << ".control";  // at restart -1 required here
          filenameout = name.str();
          filenameout = filenameout.substr(0, filenameout.length() - 8);
          break;
        }
      }

      // delete -0 at end of filename if present
      pos = filenameout.rfind('-');
      number = atoi(filenameout.substr(pos + 1).c_str());
      if (number == 0) filenameout = filenameout.substr(0, filenameout.length() - 2);
    }
  }

  return filenameout;
}
