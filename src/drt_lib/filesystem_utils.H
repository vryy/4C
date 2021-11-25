/*! \file
\level 1
\brief Util functions for the filesystem
*/

#ifndef FILESYSTEM_UTILS_H_
#define FILESYSTEM_UTILS_H_


#include <string>

namespace UTILS
{
  namespace FILESYSTEM
  {
    /*!
     * @brief Returns true if the given path is an absolute path
     *
     * @param pathname Pathname
     * @return true
     * @return false
     */
    inline bool IsAbsolutePath(const std::string& pathname) { return pathname.find('/') == 0; }

    /*!
     * @brief Returns the directory name of the filename
     *
     * @param filename filename
     * @return std::string name of the directory
     */
    inline std::string GetDirectoryName(const std::string& filename)
    {
      if (filename.find('/') == std::string::npos) return "./";

      std::size_t last_sep = filename.rfind('/');

      return filename.substr(0, last_sep + 1);
    }

    /*!
     * @brief Joins two path strings
     *
     * @param path1
     * @param path2
     * @return std::string
     */
    inline std::string Join(const std::string& path1, const std::string& path2)
    {
      // check if path2 is an absolute path
      if (path2.rfind('/', 0) == 0 or path1.length() == 0) return path2;

      std::stringstream ss;

      ss << path1;

      if (path1.back() != '/')
      {
        ss << '/';
      }

      ss << path2;

      return ss.str();
    }
  }  // namespace FILESYSTEM
}  // namespace UTILS

#endif  // FILESYSTEM_UTILS_H_