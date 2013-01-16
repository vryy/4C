/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils.H

\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef DRT_CONDITION_UTILS_H
#define DRT_CONDITION_UTILS_H


#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>


namespace LINALG
{
  class MapExtractor;
}

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
  class Node;

  namespace UTILS
  {
    // forward declaration
    class ConditionSelector;
    class DiscretizationNodeIterator;

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

    /// find all local nodes from discretization marked with condition
    /*!
      Loop all conditions of the given discretization, find the ones with the
      specified name and collect the locally owned node ids in the supplied
      set. The nodes vector is unique and ordered on output.

      \param dis : (in) Discretization
      \param condname : (in) name of condition in question
      \param nodes : (out) empty set on input, filled with nodal gids of local nodes

      \author u.kue
      \date 06/07
    */
    void FindConditionedNodes(const DRT::Discretization& dis,
                              std::string condname,
                              std::vector<int>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
                              std::string condname,
                              std::map<int, DRT::Node*>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
                                  std::string condname,
                                  std::set<int>& nodeset);

    /// find all local nodes from discretization marked with condition
    /*!
      Loop all conditions of the given discretization, find the ones with the
      specified name and collect the locally owned node ids in the suppied
      set. The nodes vector is unique and ordered on output.

      \param dis : (in) Discretization
      \param conds : (in) conditions in question
      \param nodes : (out) empty set on input, filled with nodal gids of local nodes

      \author u.kue
      \date 06/07
    */
    void FindConditionedNodes(const DRT::Discretization& dis,
                              const std::vector<DRT::Condition*>& conds,
                              std::vector<int>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
                              const std::vector<DRT::Condition*>& conds,
                              std::map<int, DRT::Node*>& nodes);

    /// find all local nodes from discretization marked with condition and
    /// put them into a vector indexed by Id of the condition
    void FindConditionedNodes(const DRT::Discretization& dis,
                              const std::vector<DRT::Condition*>& conds,
                              std::map<int, std::map<int, DRT::Node*> >& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
                              const std::vector<DRT::Condition*>& conds,
                              std::set<int>& nodeset);


    /// collect all local nodes and elements in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param elements unique map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
        std::map<int, DRT::Node*>& nodes,
        std::map<int, Teuchos::RCP<DRT::Element> >& elements,
                              const string& condname);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
                              std::map<int, DRT::Node*>& nodes,
                              std::map<int, DRT::Node*>& ghostnodes,
                              std::map< int, Teuchos::RCP<DRT::Element> >& elements,
                              std::vector<DRT::Condition*>& conds);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
                              std::map<int, DRT::Node*>& nodes,
                              std::map<int, DRT::Node*>& ghostnodes,
                              std::map< int, Teuchos::RCP<DRT::Element> >& elements,
                              const string& condname);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
                              std::map<int, std::map<int, DRT::Node*> >& nodes,
                              std::map<int, std::map<int, DRT::Node*> >& ghostnodes,
                              std::map<int, std::map< int, Teuchos::RCP<DRT::Element> > >& elements,
                              const string& condname);

    /// collect all elements in a condition including ghosts
    /*!
      \param dis discretization
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
                              std::map<int, Teuchos::RCP<DRT::Element> >& elements,
                              const string& condname);


    /// Create element map to given condition
    /*!
      \pre Geometry must be available in condition

      \param dis : (in) Discretization
      \param condname : (in) name of condition in question
      \param colmap : (in) whether we want to create a (ghosted) column map
    */
    Teuchos::RCP<Epetra_Map> ConditionElementMap(const DRT::Discretization& dis, std::string condname, bool colmap);

    /// Find all conditions with given name that all nodes of the element have in common
    /*!
      \param ele (in) the element
      \param condname (in) name of the condition to look for
      \param condition (out) all conditions that cover all element nodes
    */
    void FindElementConditions(const DRT::Element* ele, const std::string& condname, std::vector<DRT::Condition*>& condition);

    /// row map with nodes from condition
    Teuchos::RCP<Epetra_Map> ConditionNodeRowMap(const DRT::Discretization& dis,
                                                 const std::string& condname);

    /// col map with nodes from condition
    Teuchos::RCP<Epetra_Map> ConditionNodeColMap(const DRT::Discretization& dis,
                                                 const std::string& condname);

    /// generic condition node map creation
    Teuchos::RCP<Epetra_Map> ConditionMap(const DRT::Discretization& dis,
                                          const DiscretizationNodeIterator& iter,
                                          const std::string& condname);

    /// create the set of column element gids that have conditioned nodes
    /*!
      \note These are not elements from the condition geometry. Rather the
      gids of actual discretization elements are listed.
     */
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis,
                                                       const std::string& condname);


    //! create discretization from a given condition
    Teuchos::RCP<DRT::Discretization> CreateDiscretizationFromCondition(
      Teuchos::RCP<DRT::Discretization> sourcedis,    ///< discretization with condition
      const std::string&                condname,     ///< name of the condition, by which the derived discretization is identified
      const std::string&                discret_name, ///< name of the new condition
      const std::string&                element_name,  ///< name/type of the elements to be created
      const std::vector<std::string>&   conditions_to_copy  ///< list of conditions that will be copied to the new discretization
      );


    //!collects elements by labels (have to be implemented in the corresponding condition
    void CollectElementsByConditionLabel(
        const DRT::Discretization&           discret,               ///< discretization with condition
        std::map<int,std::set<int> >&        elementsByLabel,       ///< map from label to set of element GID
        const string&                        name                   ///< name of the condition
        );

  }
}

#endif
