/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the condition and material definition

\level 3


*/
/*---------------------------------------------------------------------*/


#include <string>

#include "lib_dserror.H"
#include "lib_utils_cond_and_mat_definition.H"

// Convert a string to a number, i.e. to an int or a double
// Perform the appropriate error checks
template <typename T>
T DRT::UTILS::ConvertAndValidateStringToNumber(const std::string& snumber,
    const std::string& variablename, const std::string& sectionname, int variablelength,
    bool optional)
{
  // value is set by the function stoi or stod to position of the next character in str after the
  // numerical value. Needed to check for remaining characters after string to int or double
  // conversion.
  T nnumber;
  std::size_t pos = 0;

  try
  {
    // convert to int or double, depending on type T
    nnumber = DRT::UTILS::ConvertStringToNumber<T>(snumber, pos);
  }
  catch (std::invalid_argument& e)
  {
    // in case the parameter is mandatory and no value is given
    if (!optional and snumber.empty())
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
    DRT::UTILS::ThrowErrorWrongDataType(snumber.substr(pos), nnumber, variablename, sectionname);
  }

  return nnumber;
}

template int DRT::UTILS::ConvertAndValidateStringToNumber<int>(const std::string& snumber,
    const std::string& variablename, const std::string& sectionname, int variablelength,
    bool optional);

template double DRT::UTILS::ConvertAndValidateStringToNumber<double>(const std::string& snumber,
    const std::string& variablename, const std::string& sectionname, int variablelength,
    bool optional);


// Throw an error for the wrong data type in case 'nnumber' is of type int
void DRT::UTILS::ThrowErrorWrongDataType(const std::string& snumbersubstring, int nnumber,
    const std::string& variablename, const std::string& sectionname)
{
  dserror(
      "Failed to read value '%s' while reading variable '%s' in '%s'. BACI could only read "
      "'%d', so the specified number format is probably not supported. The variable '%s' "
      "has to be an integer.",
      snumbersubstring.c_str(), variablename.c_str(), sectionname.c_str(), nnumber,
      variablename.c_str());
}

// Throw an error for the wrong data type in case 'nnumber' is of type double
void DRT::UTILS::ThrowErrorWrongDataType(const std::string& snumbersubstring, double nnumber,
    const std::string& variablename, const std::string& sectionname)
{
  dserror(
      "Failed to read value '%s' while reading variable '%s' in '%s'. BACI could only read "
      "'%f', so the specified number format is probably not supported. The variable '%s' "
      "has to be a floating point.",
      snumbersubstring.c_str(), variablename.c_str(), sectionname.c_str(), nnumber,
      variablename.c_str());
}