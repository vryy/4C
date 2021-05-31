/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the condition and material definition

\level 3


*/
/*---------------------------------------------------------------------*/


#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>

#include "drt_utils_cond_and_mat_definition.H"

// Convert a string to a number, i.e. to an int or a double
// Perform the appropriate error checks
template <typename T>
T DRT::UTILS::convertAndValidateStringToNumber(const std::string& snumber, T nnumber,
    const std::string& variablename, const std::string& sectionname, int variablelength,
    bool optional_)
{
  // value is set by the function stoi or stod to position of the next character in str after the
  // numerical value. Needed to check for remaining characters after string to int or double
  // conversion.
  std::size_t pos = 0;

  try
  {
    // convert to int or double, depending on type T
    nnumber = DRT::UTILS::convertStringToNumber<T>(snumber, pos);
  }
  catch (std::invalid_argument& e)
  {
    // in case the parameter is mandatory and no value is given
    if (!optional_ and snumber.empty())
    {
      dserror(
          "Invalid argument! No value of variable '%s' in '%s' specified. Possibly you "
          "didn't give enough input values. The variable '%s' expectes %i input values.",
          variablename.c_str(), sectionname.c_str(), variablename.c_str(), variablelength);
    }
    // any other weird input values
    else
    {
      dserror("Invalid argument! Failed to read the value '%s' of variable '%s' in '%s'.",
          snumber.c_str(), variablename.c_str(), sectionname.c_str(), variablename.c_str());
    }
  }
  // check if there are any other characters that were not converted
  if (pos != snumber.size())
  {
    DRT::UTILS::throwErrorWrongDataType<T>(snumber.substr(pos), nnumber, variablename, sectionname);
  }

  return nnumber;
}

template int DRT::UTILS::convertAndValidateStringToNumber<int>(const std::string& snumber,
    int nnumber, const std::string& variablename, const std::string& sectionname,
    int variablelength, bool optional_);

template double DRT::UTILS::convertAndValidateStringToNumber<double>(const std::string& snumber,
    double nnumber, const std::string& variablename, const std::string& sectionname,
    int variablelength, bool optional_);