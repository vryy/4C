/*----------------------------------------------------------------------*/
/*!
\file adapter_utils.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_UTILS_H
#define ADAPTER_UTILS_H

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_mapextractor.H"

namespace ADAPTER
{
  /// Helper functions to glue things together
  namespace UTILS
  {
    /// setup the extractor object to couple ndim dofs at the named condition
    /*!

      Partition the row dof map of the given discretization in two. Put this
      first ndim dofs of the nodes covered by the given condition in slot
      one. All the other dofs in slot zero.

      \author u.kue
      \date 01/08
     */
    void SetupNDimExtractor(const DRT::Discretization& dis,
                            std::string condname,
                            LINALG::MapExtractor& extractor);

    /// row map with nodes from condition
    Teuchos::RCP<Epetra_Map> ConditionNodeMap(const DRT::Discretization& dis,
                                              std::string condname);
  }
}

#endif
#endif
