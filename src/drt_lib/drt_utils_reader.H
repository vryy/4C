/*!----------------------------------------------------------------------
\file drt_utils_reader.H

\brief A collection of helper methods for readers

<pre>
\level 0

\maintainer Martin Kronbichler
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_READER_H
#define UTILS_READER_H
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#endif


namespace DRT
{
  /// collection of useful and important but less central stuff
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
    static std::string trim(const std::string& line)
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

    // remove comments, trailing and leading whitespaces
    static std::string strip_comment(const std::string& line)
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
      newline = trim(newline);

      return newline;
    }

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_READER_H
