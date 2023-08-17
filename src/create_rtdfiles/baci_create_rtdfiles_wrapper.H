/*----------------------------------------------------------------------*/
/*! \file


\brief Main routine

\author scheider

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_CREATE_RTDFILES_WRAPPER_H
#define BACI_CREATE_RTDFILES_WRAPPER_H
#include <string>

namespace DRT
{
  namespace RTD
  {
    /*!
      \brief Create a restructuredText file containing all header sections and parameters

      \param[in] string: filename for the the header parameters, should be headerreference.rst

    */
    void WriteReadTheDocsHeader(const std::string& headerdocumentationfilename);

    /*!
      \brief WCreate a restructuredText file containing all material models and parameters

      \param[in] string: filename for the the material parameters, should be materialreference.rst

    */
    void WriteReadTheDocsMaterial(const std::string& materialdocumentationfilename);

    /*!
      \brief Create a restructuredText file containing all prescribed conditions

      \param[in] string: filename for the the header parameters, should be conditionreference.rst

    */
    void WriteReadTheDocsCondition(const std::string& conditiondocumentationfilename);

    /*!
      \brief Create a restructuredText file containing all other reference sections
      at this time, it includes: RESULT DESCRIPTION, FUNCTION

      \param[in] string: filename for the the header parameters, should be furtherreference.rst

    */
    void WriteReadTheDocsVarious(const std::string& variousdocumentationfilename);


    void PrintHelpMessage();

  }  // namespace RTD
}  // namespace DRT

#endif  // CREATE_RTDFILES_WRAPPER_H
