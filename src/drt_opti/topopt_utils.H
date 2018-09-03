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


#ifndef TOPOPT_UTILS_H_
#define TOPOPT_UTILS_H_


#include <Teuchos_RCP.hpp>


// forward declarations
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
    explicit TopoptFluidCloneStrategy() {}
    /// destructor
    virtual ~TopoptFluidCloneStrategy() {}

   protected:
    /// determine element type std::string and whether element is copied or not
    bool DetermineEleType(
        DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    /// set element-specific data (material etc.)
    void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid,
        const bool isnurbs);

    /// returns conditions names to be copied (source and target name)
    std::map<std::string, std::string> ConditionsToCopy();

    /// check for correct material
    void CheckMaterialType(const int matid);

   private:
  };  // class TopoptFluidCloneStrategy

  /// prints the logo of topology optimization
  void printTopOptLogo();

  /// expands a filename
  /*
   * the given filename is expanded so that
   *    o the folder remains the same
   *    o the filename within the folder starts with the expansion
   *    o optionally the restart number is lowered by 1 (-1 vanishes)
   */
  std::string modifyFilename(const std::string& filename, const std::string& expansion,
      bool handleRestart = true, bool firstcall = false);
}  // namespace TOPOPT

#endif /* TOPOPT_UTILS_H_ */
