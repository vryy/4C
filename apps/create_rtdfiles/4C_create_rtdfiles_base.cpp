/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\level 0
*/
/*----------------------------------------------------------------------*/

#
#include "4C_config.hpp"

#include "4C_create_rtdfiles_wrapper.hpp"

#include <cstring>
#include <filesystem>
#include <iostream>


int main(int argc, char *argv[])
{
  using namespace FourC;

  MPI_Init(&argc, &argv);

  printf(
      "\n"
      "**********************************************\n"
      "*                    4C                      *\n"
      "*          Reference files creator           *\n"
      "*              for ReadTheDocs               *\n"
      "**********************************************\n\n");

  if ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)))
  {
    printf("\n\n");
    RTD::PrintHelpMessage();
    printf("\n\n");
  }
  else
  {
    std::string reference_path = (argc == 2) ? argv[1] : "reference_docs";
    if (not std::filesystem::exists(reference_path))
      std::filesystem::create_directory(reference_path);

    RTD::write_cell_type_information(reference_path + "/elementinformation.yaml");
    std::cout << "Writing cell type information to yaml finished\n";
    RTD::write_read_the_docs_header(reference_path + "/headerreference.rst");
    std::cout << "Writing headerreference.rst finished\n";
    RTD::write_read_the_docs_celltypes(reference_path + "/celltypereference.rst");
    std::cout << "Writing celltypes.rst finished\n";
    RTD::write_read_the_docs_material(reference_path + "/materialreference.rst");
    std::cout << "Writing materialreference.rst finished\n";
    RTD::write_read_the_docs_condition(reference_path + "/conditionreference.rst");
    std::cout << "Writing conditionreference.rst finished\n";
    RTD::write_read_the_docs_various(reference_path + "/furtherreference.rst");
    std::cout << "Writing furtherreference.rst finished\n";
    const std::string contactconstitivedocumentationfilename =
        reference_path + "/contactconstitutivereference.rst";
    const std::string elementdocumentationfilename = reference_path + "/elementreference.rst";
    /*
     * TODO: Other files can be written as readthedocs Reference files (e.g. Element Reference)
     */
  }

  MPI_Finalize();

  return (0);
}
