/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for readers

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef LIB_UTILS_READER_H
#define LIB_UTILS_READER_H
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <cctype>
#include <cstring>
#include <boost/algorithm/string/split.hpp>
#endif


namespace DRT
{
  namespace UTILS
  {
    /*----------------------------------------------------------------------*/
    /*!
      \brief remove all leading and trailing whitespaces from a string.

      Note: consecutive whitespaces inside the std::string will be reduced to a
      single space.

      \author u.kue
      \date 03/07
    */
    /*----------------------------------------------------------------------*/
    std::string Trim(const std::string& line);


    /**
     * Split the @p input string into multiple substrings between the @p delimiter.
     */
    std::vector<std::string> Split(const std::string& input, const std::string& delimiter);

    /*----------------------------------------------------------------------*/
    /*!
      remove comments, trailing and leading whitespaces

      @param[in] line arbitrary string
      @result same string stripped
    */
    /*----------------------------------------------------------------------*/
    std::string StripComment(const std::string& line);

    /*----------------------------------------------------------------------*/
    /*!
      convert all characters in a string into lower case
      (wrapper to the corresponding routine in boost::algorithm)

      @param[in] line arbitrary string
      @result same string in all lower case
    */
    /*----------------------------------------------------------------------*/
    std::string ToLower(const std::string& line);
  }  // namespace UTILS
}  // namespace DRT


#endif  // LIB_UTILS_READER_H
