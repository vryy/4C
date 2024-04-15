/*----------------------------------------------------------------------------*/
/*! \file
    \brief virtual DG Element interface for HDG discretization method

    \level 3
 */
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DG_ELEMENT_HPP
#define FOUR_C_LIB_DG_ELEMENT_HPP

#include "baci_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class DgElement
  {
   public:
    /*!
     \brief Destructor
    */
    virtual ~DgElement() = default;
    virtual int NumDofPerNodeAuxiliary() const = 0;

    virtual int NumDofPerElementAuxiliary() const = 0;
  };
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
