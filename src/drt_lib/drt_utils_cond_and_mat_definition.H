/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the condition and material definition

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef UTILS_COND_AND_MAT_DEFINITION_H
#define UTILS_COND_AND_MAT_DEFINITION_H


#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>

namespace DRT
{
  namespace UTILS
  {
    /*! brief Convert a string to an integer of double depending on template type T and perform
     * appropriate error checks
     *
     *  Assume a variable 'variablename' in section 'sectionname' with an value 'nnumber' of type
     * double or int (T) Convert the string snumber, extracted from the input line, into a number of
     * type T. Check for special cases, i.e. parameter is optional and assign the appropriate
     * values. Check for remaining characters that could not be extracted and wrong number formats.
     * If this is the case, throw an error message.
     */
    template <typename T>
    T convertAndValidateStringToNumber(
        const std::string&
            snumber,  ///< string read from the line to be converted to an int or double
        T nnumber,    ///< int or double variable the string shall be parsed into
        const std::string& variablename,  ///< name of the variable
        const std::string&
            sectionname,     ///< name of the section the variable and its value is read from
        int variablelength,  ///< length of the varible
        bool optional_       ///< optional flag
    );

    /*! brief  Throw an error for the wrong data type in case 'nnumber' is of type T but some other
     * type /neq T was given
     *
     *  Assume a variable 'variablename' in section 'sectionname' with an value 'nnumber' and the
     * 'snumbersubstring' that could not be converted to 'nnumber'.
     */
    template <typename T>
    void throwErrorWrongDataType(
        const std::string&
            snumbersubstring,  ///< remaining characters that could not be converted to a number
        T nnumber,             ///< number the error refers to
        const std::string& variablename,  ///< name of the variable
        const std::string&
            sectionname  ///< name of the section the variable and its value is read from
    )
    {
      dserror("You shouldn't land here - data type can be double or int.");
    }

    /// specialization for int
    template <>
    inline void throwErrorWrongDataType<int>(const std::string& snumbersubstring, int nnumber,
        const std::string& variablename, const std::string& sectionname)
    {
      dserror(
          "Failed to read value '%s' while reading variable '%s' in '%s'. BACI could only read "
          "'%d', so the specified number format is probably not supported. The variable '%s' "
          "has to be an integer.",
          snumbersubstring.c_str(), variablename.c_str(), sectionname.c_str(), nnumber,
          variablename.c_str());
    }

    /// specialization for double
    template <>
    inline void throwErrorWrongDataType<double>(const std::string& snumbersubstring, double nnumber,
        const std::string& variablename, const std::string& sectionname)
    {
      dserror(
          "Failed to read value '%s' while reading variable '%s' in '%s'. BACI could only read "
          "'%f', so the specified number format is probably not supported. The variable '%s' "
          "has to be a floating point.",
          snumbersubstring.c_str(), variablename.c_str(), sectionname.c_str(), nnumber,
          variablename.c_str());
    }

    /*! brief  Convert the string 'snumber' to number of type T and return that value number via
     * stoi or stod
     *
     *  Depending on whether T is type int or double, use stoi or stod to convert 'snumber'.
     */
    template <typename T>
    T convertStringToNumber(
        const std::string&
            snumber,       ///< string read from the line to be converted to an int or double
        std::size_t& pos)  ///< position of the next character in string after the numerical value
    {
      dserror("You shouldn't land here - data type can be double or int.");
    }

    /// specialization for int
    template <>
    inline int convertStringToNumber<int>(const std::string& snumber, std::size_t& pos)
    {
      return stoi(snumber, &pos);
    }

    /// specialization for double
    template <>
    inline double convertStringToNumber<double>(const std::string& snumber, std::size_t& pos)
    {
      return stod(snumber, &pos);
    }

  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifndef  UTILS_COND_AND_MAT_DEFINITION_H
