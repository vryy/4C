/*----------------------------------------------------------------------*/
/*!
\file ale_utils_clonestrategy.H

\brief mesh clone strategy for ale problems

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ALE_UTILS_CLONESTRATEGY_H
#define ALE_UTILS_CLONESTRATEGY_H

#include <Teuchos_RCP.hpp>
#include <vector>
#include <map>
#include <string>

// forward declarations
namespace DRT
{
 class Element;
}

// the actual classes
namespace ALE
{
  namespace UTILS
  {
    /*!
    \brief Implementation of special clone strategy for automatic generation
           of ale from a given fluid discretization

    */
    class AleCloneStrategy
    {
    public:

    /// constructor
    explicit AleCloneStrategy(){}
    /// destructor
    virtual ~AleCloneStrategy(){}

    protected:

    /// determine element type string and whether element is copied or not
    bool DetermineEleType(
        DRT::Element* actele,               ///< current element
        const bool ismyele,                 ///< true if element belongs to my proc
        std::vector<std::string>& eletype); ///< element type

    /// set element-specific data (material etc.)
    void SetElementData(
        Teuchos::RCP<DRT::Element> newele,
        DRT::Element* oldele,
        const int matid,
        const bool nurbsdis);

    /// returns conditions names to be copied (source and target name)
    std::map<std::string,std::string> ConditionsToCopy();

    /// check for correct material
    void CheckMaterialType(const int matid);

    private:

    }; // class AleCloneStrategy
  } // namespace UTILS
} // namespace ALE

#endif
