/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\author scheider

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CREATE_RTDFILES_WRAPPER_HPP
#define FOUR_C_CREATE_RTDFILES_WRAPPER_HPP
#include "4C_config.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace RTD
{
  /*!
    \brief Create a yaml file containing the cell type information

    \param[in] string: filename for the the yaml file, should be elementinformation.yaml
  */
  void write_cell_type_information(const std::string& elementinformationfilename);

  /*!
    \brief Create a restructuredText file containing all header sections and parameters

    \param[in] string: filename for the the header parameters, should be headerreference.rst
  */
  void write_read_the_docs_header(const std::string& headerdocumentationfilename);

  /*!
    \brief Create a restructuredText file containing the cell type section
      including node numbers, etc.
      TODO face IDs, etc.

    \param[in] string: filename for the the header parameters, should be celltypereference.rst

  */
  void write_read_the_docs_celltypes(const std::string& celltypedocumentationfilename);

  /*!
    \brief WCreate a restructuredText file containing all material models and parameters

    \param[in] string: filename for the the material parameters, should be materialreference.rst

  */
  void write_read_the_docs_material(const std::string& materialdocumentationfilename);

  /*!
    \brief Create a restructuredText file containing all prescribed conditions

    \param[in] string: filename for the the header parameters, should be conditionreference.rst

  */
  void write_read_the_docs_condition(const std::string& conditiondocumentationfilename);

  /*!
    \brief Create a restructuredText file containing all other reference sections
    at this time, it includes: RESULT DESCRIPTION, FUNCTION

    \param[in] string: filename for the the header parameters, should be furtherreference.rst

  */
  void write_read_the_docs_various(const std::string& variousdocumentationfilename);


  void PrintHelpMessage();

}  // namespace RTD

FOUR_C_NAMESPACE_CLOSE

#endif
