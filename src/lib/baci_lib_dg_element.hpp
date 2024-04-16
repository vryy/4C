/*----------------------------------------------------------------------------*/
/*! \file
    \brief virtual DG Element interface for HDG discretization method

    \level 3
 */
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DG_ELEMENT_HPP
#define FOUR_C_LIB_DG_ELEMENT_HPP

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

#endif
