/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_factory.H

\brief factory for scatra boundary evaluation

\level 2

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
 */
/*----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_BOUNDARY_FACTORY_H
#define SCATRA_ELE_BOUNDARY_FACTORY_H

#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_scatra.H"

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
      virtual ~ScaTraBoundaryFactory() { return; };

      //! ProvideImpl
      static ScaTraBoundaryInterface* ProvideImpl(const DRT::Element* ele,
          const enum INPAR::SCATRA::ImplType impltype, const int numdofpernode, const int numscal,
          const std::string& disname);

     private:
      //! return instance of element evaluation class depending on implementation type
      template <DRT::Element::DiscretizationType distype>
      static ScaTraBoundaryInterface* DefineProblemType(const enum INPAR::SCATRA::ImplType impltype,
          const int numdofpernode, const int numscal, const std::string& disname);
    };  // class ScaTraBoundaryFactory
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
