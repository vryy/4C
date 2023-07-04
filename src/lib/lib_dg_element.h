/*----------------------------------------------------------------------------*/
/*! \file
    \brief virtual DG Element interface for HDG discretization method

    \level 3
 */
/*----------------------------------------------------------------------------*/

#ifndef LIB_DG_ELEMENT_H
#define LIB_DG_ELEMENT_H

namespace DRT
{
  class DG_Element
  {
   public:
    /*!
     \brief Destructor
    */
    virtual ~DG_Element() {}

    virtual int NumDofPerNodeAuxiliary() const = 0;

    virtual int NumDofPerElementAuxiliary() const = 0;
  };
}  // namespace DRT



#endif  // BACI_LIB_DG_ELEMENT_H
