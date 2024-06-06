/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\author scheider

\level 0
*/

#include "4C_config_revision.hpp"

#include "4C_create_rtdfiles_wrapper.hpp"

#include "4C_comm_utils.hpp"
#include "4C_create_rtdfiles_utils.hpp"
#include "4C_discretization_fem_general_element_definition.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validconditions.hpp"
#include "4C_inpar_validcontactconstitutivelaw.hpp"
#include "4C_inpar_validmaterials.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

namespace RTD
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteCellTypeInformation(const std::string& elementinformationfilename)
  {
    // open ascii file for writing the cell type information
    std::ofstream elementinformationfile(elementinformationfilename.c_str());
    if (!elementinformationfile)
      FOUR_C_THROW("failed to open file: %s", elementinformationfilename.c_str());
    elementinformationfile << "# yaml file created using 4C version (git SHA1):\n";
    elementinformationfile << "# " << VersionControl::git_hash << "\n#\n";

    WriteYamlCellTypeInformation(elementinformationfile);
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteReadTheDocsHeader(const std::string& headerdocumentationfilename)
  {
    // open ascii file for writing all header parameters
    std::ofstream headerdocumentationfile(headerdocumentationfilename.c_str());
    if (!headerdocumentationfile)
      FOUR_C_THROW("failed to open file: %s", headerdocumentationfilename.c_str());
    headerdocumentationfile << "..\n   Created using 4C version (git SHA1):\n";
    headerdocumentationfile << "   " << VersionControl::git_hash << "\n\n";
    headerdocumentationfile << ".. _headerparameters:\n\n";
    headerdocumentationfile << "Header parameters\n";
    headerdocumentationfile << "=================\n\n";
    WriteHeaderReference(headerdocumentationfile, *Input::ValidParameters(), "");
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteReadTheDocsCelltypes(const std::string& celltypedocumentationfilename)
  {
    // open ascii file for writing all header parameters
    std::ofstream celltypeocumentationfile(celltypedocumentationfilename.c_str());
    if (!celltypeocumentationfile)
      FOUR_C_THROW("failed to open file: %s", celltypedocumentationfilename.c_str());
    celltypeocumentationfile << "..\n   Created using 4C version (git SHA1):\n";
    celltypeocumentationfile << "   " << VersionControl::git_hash << "\n\n";

    WriteCelltypeReference(celltypeocumentationfile);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteReadTheDocsMaterial(const std::string& materialdocumentationfilename)
  {
    //
    // open ascii file for writing all material parameters
    std::ofstream materialdocumentationfile(materialdocumentationfilename.c_str());
    if (!materialdocumentationfile)
      FOUR_C_THROW("failed to open file: %s", materialdocumentationfilename.c_str());
    materialdocumentationfile << "..\n   Created using 4C version (git SHA1):\n";
    materialdocumentationfile << "   " << VersionControl::git_hash << "\n\n";
    WriteMaterialReference(materialdocumentationfile, *Input::ValidMaterials());
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteReadTheDocsCondition(const std::string& conditiondocumentationfilename)
  {
    //
    // open ascii file for writing all constrains / conditions parameters
    std::ofstream conditiondocumentationfile(conditiondocumentationfilename.c_str());
    if (!conditiondocumentationfile)
      FOUR_C_THROW("failed to open file: %s", conditiondocumentationfilename.c_str());
    conditiondocumentationfile << "..\n   Created using 4C version (git SHA1):\n";
    conditiondocumentationfile << "   " << VersionControl::git_hash << "\n\n";
    WriteConditionsReference(conditiondocumentationfile, *Input::ValidConditions());

    WriteContactLawReference(conditiondocumentationfile, *Input::ValidContactConstitutiveLaws());
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteReadTheDocsVarious(const std::string& variousdocumentationfilename)
  {
    //
    // open ascii file for writing other (non header) parameters
    std::ofstream variousdocumentationfile(variousdocumentationfilename.c_str());
    if (!variousdocumentationfile)
      FOUR_C_THROW("failed to open file: %s", variousdocumentationfilename.c_str());
    variousdocumentationfile << "..\n   Created using 4C version (git SHA1):\n";
    variousdocumentationfile << "   " << VersionControl::git_hash << "\n\n";
    WriteVariousReference(variousdocumentationfile);
  }

  void PrintHelpMessage()
  {
    std::cout << "This program writes all necessary reference files for readthedocs\n";
    std::cout << "Usage:\n    create_rtd [pathanem]\n";
    std::cout << " Parameter:\n   pathname (str) path where the reference files are stored.\n";
    std::cout << "                   Default: reference_docs";
  }

}  // namespace RTD

FOUR_C_NAMESPACE_CLOSE
