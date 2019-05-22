/*----------------------------------------------------------------------*/
/*!

\brief a class to manage an enhanced discretization including varying number
of dofs per node on a fluid discretization for xwall

\level 2

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef DISCRET_XWALL_H
#define DISCRET_XWALL_H

#include "drt_discret_faces.H"
#include "drt_dserror.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include <string>
#include "Teuchos_ParameterList.hpp"

namespace DRT
{
  class DiscretizationXWall : public DiscretizationFaces
  {
   public:
    /*!
    \brief Standard Constructor

    \param name (in): name of this discretization
    \param comm (in): An epetra comm object associated with this discretization
    */
    DiscretizationXWall(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

    /*!
    \brief Destructor

    */
    virtual ~DiscretizationXWall() { return; };

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))

    Additional input nodal dof set: If the node contains more than one set of dofs, which can be
    evaluated, the number of the set needs to be given. Currently only the case for XFEM.

    \param dof (in)         : vector of dof gids (to be filled)
    \param nds (in)         : number of dofset
    \param nodaldofset (in) : number of nodal dofset
    \param node (in)        : the node
    \param element (in)     : the element (optionally)
    */
    virtual void Dof(std::vector<int>& dof, const Node* node, unsigned nds, unsigned nodaldofset,
        const Element* element = NULL) const
    {
      if (nds > 1) dserror("xwall discretization can only handle one dofset at the moment");

      dsassert(nds < dofsets_.size(), "undefined dof set");
      dsassert(havedof_, "no dofs assigned");

      std::vector<int> totaldof;
      dofsets_[nds]->Dof(totaldof, node, nodaldofset);

      if (element == NULL && element->Shape() == DRT::Element::hex8)
        dserror("element required for location vector of hex8 element");

      int size;
      if (element != NULL)
        size = std::min((int)totaldof.size(), element->NumDofPerNode(*node));
      else
        size = (int)totaldof.size();
      // only take the first dofs that have a meaning for all elements at this node
      for (int i = 0; i < size; i++) dof.push_back(totaldof.at(i));


      return;
    }
  };  // class DiscretizationXFEM
}  // namespace DRT

/// << operator
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationXWall& dis);


#endif  // #ifndef DISCRET_XWALL_H
