/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\level 0
*/
/*----------------------------------------------------------------------*/

#
#include <iostream>
#include <cstring>
#include "baci_create_rtdfiles_wrapper.H"
#include "baci_config_compile_settings.H"


int main(int argc, char *argv[])
{
  printf(
      "\n"
      "**********************************************\n"
      "*                   BACI                     *\n"
      "*          Reference files creator           *\n"
      "*              for ReadTheDocs               *\n"
      "**********************************************\n\n");

  if ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)))
  {
    printf("\n\n");
    DRT::RTD::PrintHelpMessage();
    printf("\n\n");
  }
  else
  {
    DRT::RTD::WriteReadTheDocsHeader("headerreference.rst");
    std::cout << "Writing headerreference.rst finished\n";
    DRT::RTD::WriteReadTheDocsMaterial("materialreference.rst");
    std::cout << "Writing materialreference.rst finished\n";
    DRT::RTD::WriteReadTheDocsCondition("conditionreference.rst");
    std::cout << "Writing conditionreference.rst finished\n";
    DRT::RTD::WriteReadTheDocsVarious("furtherreference.rst");
    std::cout << "Writing furtherreference.rst finished\n";
    const std::string contactconstitivedocumentationfilename = "contactconstitutivereference.rst";
    const std::string elementdocumentationfilename = "elementreference.rst";
    /*
     * TODO: Other files can be written as readthedocs Reference files (e.g. Element Reference)
     */
  }
  return (0);
}
