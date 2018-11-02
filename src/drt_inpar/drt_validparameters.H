/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.H

\brief Setup of the list of valid input parameters

\level 1

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_VALIDPARAMETERS_H
#define DRT_VALIDPARAMETERS_H

#include "../drt_io/io_pstream.H"

#include <string>
#include <iostream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DRT
{
  namespace INPUT
  {
    /// local wrapper to test multiple versions of "Yes", "YES", etc
    void BoolParameter(std::string const& paramName, std::string const& value,
        std::string const& docString, Teuchos::ParameterList* paramList);

    /// local wrapper for Teuchos::setIntParameter() that allows only integers
    void IntParameter(std::string const& paramName, int const value, std::string const& docString,
        Teuchos::ParameterList* paramList);

    /// local wrapper for Teuchos::setDoubleParameter() that allows only doubles
    void DoubleParameter(std::string const& paramName, double const& value,
        std::string const& docString, Teuchos::ParameterList* paramList);

    /// special implementation for a parameter being an arbitrary std::string
    void StringParameter(std::string const& paramName, std::string const& value,
        std::string const& docString, Teuchos::ParameterList* paramList);

    /// construct list with all parameters and documentation
    Teuchos::RCP<const Teuchos::ParameterList> ValidParameters();

    /// print all parameters that have a default value
    void PrintDefaultParameters(IO::Pstream& stream, const Teuchos::ParameterList& list);

    /// print flag sections of dat file with given list
    void PrintDatHeader(std::ostream& stream, const Teuchos::ParameterList& list,
        std::string parentname = "", bool color = false, bool comment = true);

    /// provide std::string to keep flag for printing equal signs
    ///
    /// This is needed for the NOX parameters whose key words
    /// and value have white spaces and thus '=' are inserted
    /// to distinguish them.
    std::string PrintEqualSign();

    /// set printing equal signs
    void SetPrintEqualSign(Teuchos::ParameterList& list, const bool& pes);

    /// Does the actual parameter list need to print equal signs?
    bool NeedToPrintEqualSign(const Teuchos::ParameterList& list);

  }  // namespace INPUT
}  // namespace DRT


/*! print list of valid parameters with documentation */
void PrintValidParameters();

/*! print help message */
void PrintHelpMessage();

/*! print flag sections of dat file with default flags */
void PrintDefaultDatHeader();


#endif
