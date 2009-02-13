/*----------------------------------------------------------------------*/
/*!
\file drt_validmaterials.H

\brief Setup of the list of valid materials for input

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef CCADISCRET

#ifndef DRT_VALIDMATERIALS_H
#define DRT_VALIDMATERIALS_H

#ifdef __cplusplus /* stuff which is c++ specific */

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>
#include <string>
#include <iostream>

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>


namespace DRT {
namespace INPUT {

  class MaterialDefinition;

  /// construct list with all materials and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > ValidMaterials();

  /// print all known material sections without contents
  void PrintEmptyMaterialDefinitions(
    std::ostream& stream,
    std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& condlist,
    bool color=true
    );

}
}

#endif

/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

  /*! print empty material sections */
  void PrintMaterialDatHeader();

#ifdef __cplusplus
}
#endif

#endif
#endif
