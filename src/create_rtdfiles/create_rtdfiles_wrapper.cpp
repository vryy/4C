/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\author scheider

\level 0
*/

#include "create_rtdfiles_wrapper.H"
#include "lib_dserror.H"
#include <iostream>
#include "lib_globalproblem.H"
#include "inpar_validparameters.H"
#include "inpar_validcontactconstitutivelaw.H"
#include "inpar_validconditions.H"
#include "inpar_validmaterials.H"
#include "lib_function.H"
#include "lib_elementdefinition.H"

#include "comm_utils.H"
#include "config_compile_settings.H"
#include "lib_utils_createdis.H"
#include "create_rtdfiles_utils.H"

#include "config_revision.H"

namespace DRT
{
  namespace RTD
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void WriteReadTheDocsHeader(const std::string& headerdocumentationfilename)
    {
      // open ascii file for writing all header parameters
      std::ofstream headerdocumentationfile(headerdocumentationfilename.c_str());
      if (!headerdocumentationfile)
        dserror("failed to open file: %s", headerdocumentationfilename.c_str());
      headerdocumentationfile << "..\n   Created using baci version (git SHA1):\n";
      headerdocumentationfile << "   " << BaciGitHash.c_str() << "\n\n";

      headerdocumentationfile << ".. _headerparameters:\n\n";
      headerdocumentationfile << "Header parameters\n";
      headerdocumentationfile << "=================\n\n";
      WriteHeaderReference(headerdocumentationfile, *DRT::INPUT::ValidParameters(), "");
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void WriteReadTheDocsMaterial(const std::string& materialdocumentationfilename)
    {
      //
      // open ascii file for writing all material parameters
      std::ofstream materialdocumentationfile(materialdocumentationfilename.c_str());
      if (!materialdocumentationfile)
        dserror("failed to open file: %s", materialdocumentationfilename.c_str());
      materialdocumentationfile << "..\n   Created using baci version (git SHA1):\n";
      materialdocumentationfile << "   " << BaciGitHash.c_str() << "\n\n";
      WriteMaterialReference(materialdocumentationfile, *DRT::INPUT::ValidMaterials());
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void WriteReadTheDocsCondition(const std::string& conditiondocumentationfilename)
    {
      //
      // open ascii file for writing all constrains / conditions parameters
      std::ofstream conditiondocumentationfile(conditiondocumentationfilename.c_str());
      if (!conditiondocumentationfile)
        dserror("failed to open file: %s", conditiondocumentationfilename.c_str());
      conditiondocumentationfile << "..\n   Created using baci version (git SHA1):\n";
      conditiondocumentationfile << "   " << BaciGitHash.c_str() << "\n\n";
      WriteConditionsReference(conditiondocumentationfile, *DRT::INPUT::ValidConditions());
    }

    void PrintHelpMessage()
    {
      std::cout << "This program writes all necessary reference files for readthedocs";
    }

  }  // namespace RTD
}  // namespace DRT