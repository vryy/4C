/*----------------------------------------------------------------------*/
/*! \file

\brief data structures for EAS related stuff in solid elements

\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_EAS_UTILS_H
#define SOLID_ELE_EAS_UTILS_H

namespace STR::ELEMENTS
{
  enum class EASType  // with meaningfull value for matrix size info
  {
    soh8_easnone,
    eastype_h8_9,
    eastype_h8_21,
    eastype_sh8_7,
    eastype_sh18_9,
    eastype_undefined

  };
}  // namespace STR::ELEMENTS

#endif  // SOLID_ELE_EAS_UTILS_H
