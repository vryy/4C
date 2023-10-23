/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of string helper methods for namespace DRT::UTILS

\level 1

*/

#include "baci_io_utils_reader.H"

namespace DRT
{
  namespace UTILS
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    std::string Trim(const std::string& line)
    {
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
      return boost::algorithm::trim_all_copy(boost::algorithm::replace_all_copy(line, "\t", " "));
#else
      std::istringstream t;
      std::string s;
      std::string newline;
      t.str(line);
      while (t >> s)
      {
        newline.append(s);
        newline.append(" ");
      }
      if (newline.size() > 0) newline.resize(newline.size() - 1);
      return newline;
#endif
    }

    std::vector<std::string> Split(const std::string& input, const std::string& delimiter)
    {
      std::vector<std::string> return_value{};
      return boost::split(return_value, input, boost::is_any_of(delimiter));
    }

    // remove comments, trailing and leading whitespaces
    std::string StripComment(const std::string& line)
    {
      std::string newline = line;

      // remove comments

      std::string::size_type loc = line.find("//");
      if (loc != std::string::npos)
      {
        newline = newline.substr(0, loc);
      }

      // remove trailing and leading whitespaces
      // compact internal whitespaces
      newline = Trim(newline);

      return newline;
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    std::string ToLower(const std::string& line) { return boost::algorithm::to_lower_copy(line); }

  }  // namespace UTILS
}  // namespace DRT
