/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the linedefinition

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef LIB_UTILS_LINEDEFINITION_H
#define LIB_UTILS_LINEDEFINITION_H

#include <memory>
#include <string>
#include <cstdlib>
#include <istream>

#include "utils_exceptions.H"


namespace DRT
{
  namespace UTILS
  {
    /*----------------------------------------------------------------------*/
    /*!
      \brief Return the data type as a string
    */
    /*----------------------------------------------------------------------*/
    template <class T>
    std::string StringFromDataType();

    /// specialisation for int
    template <>
    inline std::string StringFromDataType<int>()
    {
      return "integer";
    }

    /// specialisation for double
    template <>
    inline std::string StringFromDataType<double>()
    {
      return "double";
    }

    /// specialisation for string
    template <>
    inline std::string StringFromDataType<std::string>()
    {
      return "string";
    }

    /*----------------------------------------------------------------------*/
    /*!
      \brief Checks if any characters except comments and whitespaces were not correctly parsed from
      'stream' to a varible of type 'correctDataType'. If true, writes an error message containing
      the correct data type.
    */
    /*----------------------------------------------------------------------*/
    void CheckStreamForUnparsedCharacters(
        std::istream& stream, std::string variableName, std::string correctDataType)
    {
      if (!(stream.eof()))  // otherwise peek writes a failbit
      {
        if ((stream.peek() != ' ') && (stream.peek() != '/'))  // exclude whitespaces and comments
        {
          dserror(
              "Data type of component %s doesn't match the specified value. Needs to be of type %s",
              variableName.c_str(), correctDataType.c_str());
        }
      }
    }
  }  // namespace UTILS
}  // namespace DRT

#endif  // LIB_UTILS_LINEDEFINITION_H
