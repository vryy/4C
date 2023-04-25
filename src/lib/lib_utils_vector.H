/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef LIB_UTILS_VECTOR_H
#define LIB_UTILS_VECTOR_H

namespace DRT
{
  namespace UTILS
  {
    /*!
     * \brief sort elements and remove duplicates from vector of integers
     *
     * @param[in/out] vector            sort elements in this vector and remove duplicates
     */
    inline void SortAndRemoveDuplicateVectorElements(std::vector<int>& vector)
    {
      std::sort(vector.begin(), vector.end());
      vector.erase(std::unique(vector.begin(), vector.end()), vector.end());
    }
  }  // namespace UTILS
}  // namespace DRT


#endif  // LIB_UTILS_VECTOR_H
