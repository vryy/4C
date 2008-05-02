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

    /// std unary function version of Epetra_Map::MyGID()
    struct MyGID : public std::unary_function<int, bool>
    {
      const Epetra_Map* emap_;
      MyGID(const Epetra_Map* emap) : emap_(emap) {}
      bool operator()(int gid) const
        {
          return emap_->MyGID(gid);
        }
    };

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

    /// create ale discretization parallel to the fluid one
    void CreateAleDiscretization();

    //! create discretization from a given condition
    Teuchos::RCP<DRT::Discretization> CreateDiscretizationFromCondition(
      Teuchos::RCP<DRT::Discretization> cutterdis,    ///< condition with discretization
      const std::string&                condname,     ///< name of the condition, by which the derived discretization is identified
      const std::string&                discret_name, ///< name of the new condition
      const std::string&                element_name  ///< name/type of the elements to be created
      );
  }
}

#endif
#endif
