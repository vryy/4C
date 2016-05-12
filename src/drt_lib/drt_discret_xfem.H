/*!----------------------------------------------------------------------
\file drt_discret_xfem.H

\brief Declaration a class to manage one discretization with changing dofs

<pre>
\brief Declaration
\level 1
\maintainer Christoph Ager
            http://www.lnm.mw.tum.de
            089 - 289 - 15249
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef DISCRET_XFEM_H
#define DISCRET_XFEM_H

#include "drt_discret_faces.H"

/*!
\brief Declaration DRT: namespace of the baci discretization module

*/
namespace DRT
{

/*!
\brief Declaration A class to manage a discretization in parallel with changing dofs

*/
class DiscretizationXFEM : public DiscretizationFaces
{
public:

  /*!
  \brief Declaration Standard Constructor

  \param name (in): name of this discretization
  \param comm (in): An epetra comm object associated with this discretization
  */
  DiscretizationXFEM(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

  /*!
  \brief Declaration Destructor

  */
  virtual ~DiscretizationXFEM()
  {
    return;
  };

  /*!
  \brief Declaration Complete construction of a discretization  (Filled()==true NOT prerequisite)

  This call is done at the initial state of the discretisation, therefore the initial dofset
  is stored!

  After adding or deleting nodes or elements or redistributing them in parallel,
  or adding/deleting boundary conditions, this method has to be called to (re)construct
  pointer topologies.<br>
  It builds in this order:<br>
  - row map of nodes
  - column map of nodes
  - row map of elements
  - column map of elements
  - pointers from elements to nodes
  - pointers from nodes to elements
  - assigns degrees of freedoms
  - map of element register classes
  - calls all element register initialize methods
  - build geometries of all Dirichlet and Neumann boundary conditions

  \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                       assigning of degrees of freedoms to nodes and
                                       elements.
  \param initelements (in) : if true, build element register classes and call Initialize()
                             on each type of finite element present
  \param doboundaryconditions (in) : if true, build geometry of boundary conditions
                                     present.

  \note In order to receive a fully functional discretization, this method must be called
        with all parameters set to true (at least once). The parameters though can be
        used to turn off specific tasks to allow for more flexibility in the
        construction of a discretization, where it is known that this method will
        be called more than once.

  \note Sets Filled()=true
  */
  virtual int InitialFillComplete(bool assigndegreesoffreedom = true,
                           bool initelements = true,
                           bool doboundaryconditions = true);

  /// Export Vector with initialdofrowmap (all nodes have one dofset) - to Vector with all active dofs
  void ExportInitialtoActiveVector(
      Teuchos::RCP<const Epetra_Vector>& initialvec,
      Teuchos::RCP<Epetra_Vector>& activevec
      );


  /*!
  \brief Declaration Get the gid of all initial dofs of a node.

  Ask the initial DofSet for the gids of the dofs of this node. The
  required vector is created and filled on the fly. So better keep it
  if you need more than one dof gid.
  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  \param nds (in)       : number of dofset
  \param node (in)      : the node
  */
  virtual std::vector<int> InitialDof(unsigned nds, const Node* node) const
  {
    dsassert(nds<initialdofsets_.size(),"undefined dof set");
    dsassert(initialized_,"no initial dofs assigned");
    return initialdofsets_[nds]->Dof(node);
  }

  /*!
  \brief Declaration Get the gid of all initial dofs of a node.

  Ask the initial DofSet for the gids of the dofs of this node. The
  required vector is created and filled on the fly. So better keep it
  if you need more than one dof gid.
  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  \param node (in)      : the node
  */
  virtual std::vector<int> InitialDof(const Node* node) const
  { dsassert(initialdofsets_.size()==1,"expect just one dof set");
    dsassert(initialized_,"no initial dofs assigned");
    return InitialDof(0,node);
  }

  /*!
  \brief Declaration Get the gid of all initial dofs of a node.

  Ask the initial DofSet for the gids of the dofs of this node. The
  required vector is created and filled on the fly. So better keep it
  if you need more than one dof gid.
  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  \param nds (in)       : number of dofset
  \param node (in)      : the node
  \param lm (in/out)    : lm vector the dofs are appended to
  */

  virtual void InitialDof(unsigned nds, const Node* node, std::vector<int>& lm) const
  { dsassert(nds<initialdofsets_.size(),"undefined dof set");
    dsassert(initialized_,"no initial dofs assigned");
    initialdofsets_[nds]->Dof(node,lm); }

  /*!
  \brief Declaration Get the gid of all initial dofs of a node.

  Ask the initial DofSet for the gids of the dofs of this node. The
  required vector is created and filled on the fly. So better keep it
  if you need more than one dof gid.
  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  \param node (in)      : the node
  \param lm (in/out)    : lm vector the dofs are appended to
  */
  virtual void InitialDof(const Node* node, std::vector<int>& lm) const
  { dsassert(initialdofsets_.size()==1,"expect just one dof set");
    dsassert(initialized_,"no initial dofs assigned");
    InitialDof((unsigned)0,node,lm); }


  /// Access to initial dofset
  const DofSet& InitialDofSet(unsigned int nds=0) const
  {
    Initialized();
    return *initialdofsets_[nds];
  }
  /*!
  \brief Declaration Get initial degree of freedom row map (Initialized()==true prerequisite)

  Return ptr to initial degree of freedom row distribution map of this discretization.
  If it does not exist yet, build it.

  - Initialized()==true prerequisite

  */
  virtual const Epetra_Map* InitialDofRowMap(unsigned nds=0) const;

  /*!
  \brief Declaration Get initial degree of freedom column map (Initialized()==true prerequisite)

  Return ptr to initial degree of freedom column distribution map of this discretization.
  If it does not exist yet, build it.

  - Initialized()==true prerequisite

  */
  virtual const Epetra_Map* InitialDofColMap(unsigned nds=0) const;

  /// checks if Discretization is initialized
  bool Initialized() const;


  /*!
  \brief Declaration Set a reference to a data vector

  Using this method, a reference to a vector can
  be supplied to the discretization. The elements can access
  this vector by using the name of that vector.
  The method expects state to be either of dof row map or of
  dof column map.
  If the vector is supplied in DofColMap() a reference to it will be stored.
  If the vector is NOT supplied in DofColMap(), but in DofRowMap(),
   a vector with column map is allocated and the supplied vector is exported to it.
  Everything is stored/referenced using Teuchos::RCP.

  \param nds (in): number of dofset
  \param name (in): Name of data
  \param state (in): vector of some data

  \note This class will not take ownership or in any way modify the solution vector.
  */
  virtual void SetInitialState(unsigned nds,const std::string& name,Teuchos::RCP<const Epetra_Vector> state);


private:

  ///Store Initial Dofs
  void StoreInitialDofs();

  /*!
  ///Extend initialdofrowmap
  \param srcmap (in) : Sourcemap used as base
  \param numdofspernode (in) : Number of degrees of freedom per node
  \param numdofsets (in) : Number of XFEM-Dofsets per node
  \param uniquenumbering (in) : Assign unique number to additional dofsets

  */
  Teuchos::RCP<Epetra_Map> ExtendMap(const Epetra_Map* srcmap,
                                     int numdofspernodedofset,
                                     int numdofsets,
                                     bool uniquenumbering);

  /// initial set of dofsets
  std::vector<Teuchos::RCP<DofSet> >             initialdofsets_;

  /// bool if discretisation is initialized
  bool initialized_;

  /// full (with all reserved dofs) dof row map of initial state
  Teuchos::RCP<Epetra_Map> initialfulldofrowmap_;

  /// permuted (with duplicated gids of first dofset - to all other dofsets) dof row map of initial state
  Teuchos::RCP<Epetra_Map> initialpermdofrowmap_;


}; // class DiscretizationXFEM
} // namespace DRT

#endif  // #ifndef DISCRET_XFEM_H
