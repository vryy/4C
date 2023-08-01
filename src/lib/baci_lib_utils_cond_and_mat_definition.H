/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the condition and material definition

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_COND_AND_MAT_DEFINITION_H
#define BACI_LIB_UTILS_COND_AND_MAT_DEFINITION_H


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
    T ConvertAndValidateStringToNumber(
        const std::string&
            snumber,  ///< string read from the line to be converted to an int or double
        const std::string& variablename,  ///< name of the variable
        const std::string&
            sectionname,     ///< name of the section the variable and its value is read from
        int variablelength,  ///< length of the varible
        bool optional_       ///< optional flag
    );


    /*! brief  Throw an error for the wrong data type in case 'nnumber' is of type int but some
     * other type /neq T was given
     *
     *  Assume a variable 'variablename' in section 'sectionname' with an value 'nnumber' and the
     * 'snumbersubstring' that could not be converted to 'nnumber'.
     */
    void ThrowErrorWrongDataType(const std::string& snumbersubstring, int nnumber,
        const std::string& variablename, const std::string& sectionname);

    /*! brief  Throw an error for the wrong data type in case 'nnumber' is of type double but some
     * other type /neq T was given
     *
     *  Assume a variable 'variablename' in section 'sectionname' with an value 'nnumber' and the
     * 'snumbersubstring' that could not be converted to 'nnumber'.
     */
    void ThrowErrorWrongDataType(const std::string& snumbersubstring, double nnumber,
        const std::string& variablename, const std::string& sectionname);


    /*! brief  Convert the string 'snumber' to number of type T and return that value number via
     * stoi or stod
     *
     *  Depending on whether T is type int or double, use stoi or stod to convert 'snumber'.
     */
    template <typename T>
    T ConvertStringToNumber(
        const std::string&
            snumber,        ///< string read from the line to be converted to an int or double
        std::size_t& pos);  ///< position of the next character in string after the numerical value

    /// specialization for int
    template <>
    inline int ConvertStringToNumber<int>(const std::string& snumber, std::size_t& pos)
    {
      return stoi(snumber, &pos);
    }

    /// specialization for double
    template <>
    inline double ConvertStringToNumber<double>(const std::string& snumber, std::size_t& pos)
    {
      return stod(snumber, &pos);
    }

  }  // namespace UTILS
}  // namespace DRT

#endif  // LIB_UTILS_COND_AND_MAT_DEFINITION_H
