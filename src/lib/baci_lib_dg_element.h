/*----------------------------------------------------------------------------*/
/*! \file
    \brief virtual DG Element interface for HDG discretization method

    \level 3
 */
/*----------------------------------------------------------------------------*/

#ifndef BACI_LIB_DG_ELEMENT_H
#define BACI_LIB_DG_ELEMENT_H

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class DG_Element
  {
   public:
    /*!
     \brief Destructor
    */
    virtual ~DG_Element() = default;
    virtual int NumDofPerNodeAuxiliary() const = 0;

    virtual int NumDofPerElementAuxiliary() const = 0;
  };
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // BACI_LIB_DG_ELEMENT_H
