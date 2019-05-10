/*----------------------------------------------------------------------*/
/*!
\file sti_clonestrategy.H

\brief strategy for cloning thermo discretization from scatra discretization

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_CLONESTRATEGY_H
#define STI_CLONESTRATEGY_H

#include <Teuchos_RCP.hpp>

// forward declaration
namespace DRT
{
  class Element;
}

namespace STI
{
  /*!
  \brief strategy for cloning thermo discretization from scatra discretization

  For a scatra-thermo interaction problem, the thermo discretization is obtained through cloning
  from the scatra discretization. For convenience, we solve the thermo field as a second scatra
  field, hence the thermo discretization is in fact another scatra discretization. Unlike in other
  instances where a scatra discretization is obtained through cloning from another discretization,
  the physical implementation type of the cloned elements is set directly within the strategy and
  not in a subsequent loop after the cloning itself. We can do this because the physical
  implementation type of the cloned elements is unique and known a priori.

  \author fang (fang@lnm.mw.tum.de)
  \date 04/15
  */

  class ScatraThermoCloneStrategy
  {
   public:
    //! constructor
    explicit ScatraThermoCloneStrategy() { return; };

    //! destructor
    virtual ~ScatraThermoCloneStrategy() { return; };

    //! check material of cloned element
    void CheckMaterialType(const int matid  //! material of cloned element
    );

    //! return map with original names of conditions to be cloned as key values, and final names of
    //! cloned conditions as mapped values
    std::map<std::string, std::string> ConditionsToCopy();

    //! decide whether element should be cloned or not, and if so, determine type of cloned element
    bool DetermineEleType(DRT::Element* actele,  //! current element on source discretization
        const bool ismyele,                      //! ownership flag
        std::vector<std::string>& eletype        //! vector storing types of cloned elements
    );

    //! provide cloned element with element specific data
    void SetElementData(
        Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
        DRT::Element* oldele,               //! current element on source discretization
        const int matid,                    //! material of cloned element
        const bool isnurbs                  //! nurbs flag
    );
  };  // class ScatraThermoCloneStrategy
}  // namespace STI
#endif
