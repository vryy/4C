/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.H

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
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


/// Namespace for Baci enums of input parameters
///
/// This namespace should one day contain all the enums now still
/// defined in the various old C-headers (enums.h, dynamic.h,...)
/// The enums are used to represent the different values of std::string
/// input variables within Baci. The idea behind this new namespace
/// is twofold: first, of course, the old C-headers should go away
/// one day and we have to establish an alternative for them.
/// Second, it makes sense to define the valid parameters on a
/// global level and not within the different types of problems.
/// Of course, then drt_validparameters.H has to be included
/// wherever we want to work with the enums!
/// (a.popp 04/08)
/// \date 11/08
namespace INPAR
{
}



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
