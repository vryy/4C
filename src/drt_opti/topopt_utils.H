/*!------------------------------------------------------------------------------------------------*
\file topopt_utils.H

\brief collection of functions in namespace topopt

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_UTILS_H_
#define TOPOPT_UTILS_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;
using namespace Teuchos;

namespace DRT
{
  class Element;
}

namespace TOPOPT
{
  /*!
  \brief implementation of special clone strategy for automatic generation
  of optimization field from a given fluid discretization

  This class is nearly unchanged adopted from the according scatra class. It
  is called from the discretization and should actually be done otherwise
  (see the comment there)

  \author winklmaier
  \date 12/11
  */
  class TopoptFluidCloneStrategy
  {
  public:

    /// constructor
    explicit TopoptFluidCloneStrategy(){}
    /// destructor
    virtual ~TopoptFluidCloneStrategy() {}

  protected:

    /// determine element type string and whether element is copied or not
    bool DetermineEleType(DRT::Element* actele, const bool ismyele, std::vector<string>& eletype);

    /// set element-specific data (material etc.)
    void SetElementData(
        Teuchos::RCP<DRT::Element> newele,
        DRT::Element* oldele,
        const int matid,
        const bool isnurbs);

    /// returns conditions names to be copied (source and target name)
    map<std::string,std::string> ConditionsToCopy();

    /// check for correct material
    void CheckMaterialType(const int matid);

  private:

  }; // class ScatraFluidCloneStrategy
  //TODO check this class! It is just a copy of the scatra equivalent (winklmaier)

  // prints the logo of topology optimization
  void printTopOptLogo();


}

#endif /* TOPOPT_UTILS_H_ */
#endif  // #ifdef CCADISCRET
