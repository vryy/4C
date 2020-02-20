/*----------------------------------------------------------------------*/
/*! \file
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

    /*!
    \brief Special implementation for a parameter being an arbitrary std::string

    The method Teuchos::setNumericStringParameter() cannot be used for arbitrary
    std::string parameters, since the validate() method of the underlying
    AnyNumberParameterEntryValidator always tries to convert a given std::string to DOUBLE(s)!
    This may cause error messages in valgrind.
    Thus, for arbitrary std::strings, such as needed for specifying a file or solver name, for
    instance, this method which uses a StringValidator has to be used!

    @param[in] paramName Name of parameter to be put into the parameter list
    @param[in] value Value of the parameter
    @param[in] docString Documentation of the parameter
    @param[in/out] paramList Parameter list (to be filled with <paramName,Value,docString>)
    */
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
