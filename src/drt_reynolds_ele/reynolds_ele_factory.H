/*--------------------------------------------------------------------------*/
/*!
\file reynolds_ele_factory.H

\brief Factory of Reynolds elements

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SRC_DRT_REYNOLDS_ELE_REYNOLDS_ELE_FACTORY_H_
#define SRC_DRT_REYNOLDS_ELE_REYNOLDS_ELE_FACTORY_H_

#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace ELEMENTS
{

// forward declaration
class ReynoldsEleInterface;

class ReynoldsFactory
{
public:
  //! ctor
  ReynoldsFactory(){ return; }

  //! dtor
  virtual ~ReynoldsFactory(){ return; }

  //! ProvideImpl
  static ReynoldsEleInterface* ProvideImpl(
    DRT::Element::DiscretizationType distype,
    const std::string& disname);

private:
  //! define ReynoldsEle instances dependent on problem
  template<DRT::Element::DiscretizationType distype,int probdim>
  static ReynoldsEleInterface* DefineProblemType(
    const std::string& disname);


}; // end class ReynoldsFactory

} // end ELEMENTS

} // end DRT


#endif /* SRC_DRT_REYNOLDS_ELE_REYNOLDS_ELE_FACTORY_H_ */
