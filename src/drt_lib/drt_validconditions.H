/*----------------------------------------------------------------------*/
/*!
\file drt_validconditions.H

\brief Setup of the list of valid conditions for input

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef DRT_VALIDCONDITIONS_H
#define DRT_VALIDCONDITIONS_H

#ifdef __cplusplus /* stuff which is c++ specific */

#include <vector>
#include <string>
#include <iostream>

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>


namespace DRT {
namespace INPUT {

  class ConditionDefinition;

  /// construct list with all conditions and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > ValidConditions();

  /// print all known condition sections without contents
  void PrintEmptyConditionDefinitions(std::ostream& stream,
                                      std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

}
}

#endif

/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

  /*! print empty condition sections */
  void PrintConditionDatHeader();

#ifdef __cplusplus
}
#endif

#endif
#endif
