/*---------------------------------------------------------------------*/
/*! \file

\brief Contains ONLY lists of enumerators and is supposed to be included
       into the header files, if necessary and till the C++11 standard is
       available in BACI.


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_MATRIX_VECTOR_ENUMS_H
#define BACI_LIB_UTILS_MATRIX_VECTOR_ENUMS_H

#include <string>

namespace DRT::UTILS
{
  /// enum, which specifies the desired global vector initialization during creation
  enum class VecInitType
  {
    zero,               ///< fill the vector with zeros
    last_time_step,     ///< use the last converged time step state
    init_current_state  ///< use the current state
  };

}  // namespace DRT::UTILS



#endif  // LIB_UTILS_MATRIX_VECTOR_ENUMS_H
