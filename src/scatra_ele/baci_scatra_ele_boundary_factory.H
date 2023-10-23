/*----------------------------------------------------------------------*/
/*! \file

\brief factory for scatra boundary evaluation

\level 2

 */
/*----------------------------------------------------------------------*/


#ifndef BACI_SCATRA_ELE_BOUNDARY_FACTORY_H
#define BACI_SCATRA_ELE_BOUNDARY_FACTORY_H

#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraBoundaryInterface;

    class ScaTraBoundaryFactory
    {
     public:
      //! ctor
      ScaTraBoundaryFactory() { return; };

      //! dtor
      virtual ~ScaTraBoundaryFactory() = default;

      //! ProvideImpl
      static ScaTraBoundaryInterface* ProvideImpl(const DRT::Element* ele,
          const enum INPAR::SCATRA::ImplType impltype, const int numdofpernode, const int numscal,
          const std::string& disname);

     private:
      //! return instance of element evaluation class depending on implementation type
      template <DRT::Element::DiscretizationType distype, int probdim>
      static ScaTraBoundaryInterface* DefineProblemType(const enum INPAR::SCATRA::ImplType impltype,
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraBoundaryFactory
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
