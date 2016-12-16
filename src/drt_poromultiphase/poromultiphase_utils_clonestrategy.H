/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_utils_clonestrategy.H

 \brief utils methods for cloning the porofluid discretization


   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_CLONESTRATEGY_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_CLONESTRATEGY_H_

#include <Teuchos_RCP.hpp>

namespace DRT
{
class Element;
}

namespace POROMULTIPHASE
{
namespace UTILS
{
/*!
\brief implementation of special clone strategy for automatic generation
       of scatra from a given fluid discretization

 */
class PoroFluidMultiPhaseCloneStrategy
{
public:

  /// constructor
  explicit PoroFluidMultiPhaseCloneStrategy(){}
  /// destructor
  virtual ~PoroFluidMultiPhaseCloneStrategy() {}

  /// returns conditions names to be copied (source and target name)
  virtual std::map<std::string,std::string> ConditionsToCopy();

protected:

  /// determine element type std::string and whether element is copied or not
  virtual bool DetermineEleType(
      DRT::Element* actele,
      const bool ismyele,
      std::vector<std::string>& eletype);

  /// set element-specific data (material etc.)
  void SetElementData(
      Teuchos::RCP<DRT::Element> newele,
      DRT::Element* oldele,
      const int matid,
      const bool isnurbs);

  /// check for correct material
  void CheckMaterialType(const int matid);

private:

}; // class PoroFluidMultiPhaseCloneStrategy

} // namespace UITLS
} // namespace POROMULTIPHASE

#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_CLONESTRATEGY_H_ */
