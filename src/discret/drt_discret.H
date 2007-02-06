/*!----------------------------------------------------------------------
\file drt_discret.H

\class DRT::Discretization

\brief a class to manage one discretization

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DISCRET_H
#define DISCRET_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
using namespace std;
using namespace Teuchos;

#include "drt_element.H"
#include "drt_elementregister.H"
#include "drt_node.H"
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class Discretization

\brief A class to manage a discretization in parallel

The \ref DRT::Discretization class supports the ostream& operator <<

\author gee (gee@lnm.mw.tum.de)
*/
class Discretization
{
public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param comm : An epetra comm object associated with this discretization
  */
  Discretization(const string name, RefCountPtr<Epetra_Comm> comm);


  /*!
  \brief Destructor

  */
  virtual ~Discretization();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Get communicator associated with this class
  */
  virtual inline const Epetra_Comm& Comm() const {return *comm_;}

  /*!
  \brief Get flag indicating whether FillComplete() has been called
  */
  virtual inline bool Filled() const {return filled_;}

  /*!
  \brief Get name of this discretization
  */
  virtual inline const string& Name() const {return name_;}

  /*!
  \brief Get flag indicating whether degrees of freedom where assigned
  
  Degrees of freedom need to be assigned using AssignDegreesOfFreedom()
  before any calculations using this discretization can be made
  */
  virtual inline bool HaveDofs() const {return havedof_;}

  /*!
  \brief Get degree of freedom row map (Filled()==true prerequisite)
  
  Return ptr to degree of freedom row distribution map of this discretization.
  If it does not exist yet, build it.
  
  \note - Filled()==true prerequisite
        - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  
  */
  virtual const Epetra_Map* DofRowMap();

  /*!
  \brief Get degree of freedom column map (Filled()==true prerequisite)
  
  Return ptr to degree of freedom column distribution map of this discretization.
  If it does not exist yet, build it.
  
  \note - Filled()==true prerequisite
        - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
  
  */
  virtual const Epetra_Map* DofColMap();

  /*!
  \brief Get map associated with the distribution of the ownership of nodes
         (Filled()==true prerequisite)

  This map includes all nodes stored on this proc and also owned by this proc.
  This map is non-ambiguous, meaning that it is a non-overlapping map.

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeRowMap() const;

  /*!
  \brief Get map associated with the distribution of nodes including ghosted nodes
         (Filled()==true prerequisite)

  This map includes all nodes stored on this proc including any ghosted nodes
  This map is ambiguous, meaning that it is an overlapping map

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeColMap() const;
  /*!
  \brief Get map associated with the distribution of the ownership of elements
         (Filled()==true prerequisite)
  
  This map includes all elements stored on this proc and also owned by this proc.
  This map is non-ambiguous, meaning that it is a non-overlapping map.
  
  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* ElementRowMap() const;
  /*!
  \brief Get map associated with the distribution of elements including ghosted elements
         (Filled()==true prerequisite)
  
  This map includes all elements stored on this proc including any ghosted elements
  This map is ambiguous, meaning that it is an overlapping map
  
  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* ElementColMap() const;

  /*!
  \brief Get global number of elements (true number of total elements)
         (Filled()==true prerequisite)
  
  This is a collective call
  */
  virtual int NumGlobalElements() const; 
  
  /*!
  \brief Get processor local number of elements owned by this processor
         (Filled()==true prerequisite)
  */
  virtual int NumMyRowElements() const;

  /*!
  \brief Get processor local number of elements including ghost elements
         (Filled()==true NOT prerequisite)
  */
  virtual int NumMyColElements() const;
  
  /*!
  \brief Get global number of nodes (true number of total nodes without ghosting)
         (Filled()==true prerequisite)
  */
  virtual int NumGlobalNodes() const;

  /*!
  \brief Get processor local number of nodes owned by this processor
         (Filled()==true prerequisite)
  */
  virtual int NumMyRowNodes() const;
  
  /*! 
  \brief Get processor local number of nodes including ghost nodes
         (Filled()==true NOT prerequisite)
  */
  virtual int NumMyColNodes() const;
  
  /*!
  \brief Query whether an Element with global id gid is stored on this proc
  
  */
  virtual bool HaveGlobalElement(int gid) const;

  /*!
  \brief Get the element with global id gid (Filled()==true NOT prerequisite)
  
  Returns the element with global row id gid if element is on this proc.
  Will return row or column element, ghosted or not.
  This is an individual call
  
  \return Adress of element if element is owned by calling proc, returns NULL
          otherwise
  */
  virtual DRT::Element* gElement(int gid) const;

  /*!
  \brief Get the element with local row id lid (Filled()==true prerequisite)
  
  Returns the element with local row index lid.
  Will not return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of element if element is owned by calling proc
  */
  virtual inline DRT::Element* lRowElement(int lid) const
  { 
    #ifdef DEBUG
    if (!Filled()) dserror("DRT::Discretization::lRowElement: Filled() != true");
    #endif
    return elerowptr_[lid];
  }

  /*!
  \brief Get the element with local column id lid (Filled()==true prerequisite)
  
  Returns the element with local column index lid.
  Will also return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of element if element is stored by calling proc
  */
  virtual inline DRT::Element* lColElement(int lid) const
  {
  #ifdef DEBUG
    if (!Filled()) dserror("DRT::Discretization::lColElement: Filled() != true");
  #endif
    return elecolptr_[lid];
  }

  /*!
  \brief Query whether a Node with global id gid is stored on this proc
  
  */
  virtual bool HaveGlobalNode(int gid) const;

  /*!
  \brief Get the node with global row id gid (Filled()==true NOT prerequisite)
  
  Returns the node with global row id gid if node is on this proc.
  Will return row or column node, ghosted or  not.
  This is an individual call
  
  \return Adress of node if node is stored on calling proc
  */
  virtual DRT::Node* gNode(int gid) const;

  /*!
  \brief Get the node with local row id lid (Filled()==true prerequisite)
  
  Returns the node with local row index lid.
  Will not return any ghosted node.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of node if node is owned and stored by calling proc
  */
  virtual inline DRT::Node* lRowNode(int lid) const
  {
  #ifdef DEBUG
    if (!Filled()) dserror("DRT::Discretization::lRowNode: Filled() != true");
  #endif
    return noderowptr_[lid];
  }

  /*!
  \brief Get the node with local column id lid (Filled()==true prerequisite)
  
  Returns the node with local column index lid.
  Will return any node stored on this proc.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of node if node is stored by calling proc
  */
  virtual inline DRT::Node* lColNode(int lid) const
  {
  #ifdef DEBUG
    if (!Filled()) dserror("DRT::Discretization::lColNode: Filled() != true");
  #endif
    return nodecolptr_[lid];
  }

  /*!
  \brief Print this discretization to stdout (Filled()==true NOT prerequisite)
  
  \note This is a collective call
  */
  virtual void Print(ostream& os) const;

  //@}

  // @{ \name Construction methods

  /*!
  \brief Set an Epetra_Comm object which can be either 
         Epetra_SerialComm or Epetra_MpiComm (Filled()==true NOT prerequisite)
  
  \note Sets Filled()=false
  */
  virtual void SetComm(RefCountPtr<Epetra_Comm> comm)
  { filled_=false; comm_ = comm; return; }

  /*!
  \brief Add an element to the discretization (Filled()==true NOT prerequisite)
  
  The discretization takes ownership of the added element.
  Note that if an element with the same Id() exists, it will be 
  deleted and replaced by the new one.

  \note Sets Filled()=false
  */
  virtual void AddElement(RefCountPtr<DRT::Element> ele);

  /*!
  \brief Add a node to the discretization  (Filled()==true NOT prerequisite)
  
  The discretization takes ownership of the added node.
  Note that if a node with the same Id() exists, it will be 
  deleted and replaced by the new one.

  \note Sets Filled()=false
  */
  virtual void AddNode(RefCountPtr<DRT::Node> node);

  /*!
  \brief Complete construction of this discretization  (Filled()==true NOT prerequisite)
  
  after adding or deleting nodes or elements or redistributing them in parallel, 
  this method has to be called to (re)construct pointer topologies
  
  \note Sets Filled()=true
  */
  virtual int FillComplete();

  /*!
  \brief Set a condition with a certain name (Filled()==false on exit)
  
  Store a condition with a certain name in the discretization. The name need not
  be unique, meaning multiple conditions with the same name can be stored.
  Conditions can then be accessed with the GetCondition methods.
  
  \note Conditions attached to the discretization have to be 
        completely redundant meaning that nodal cloud in the
        condition is the same on each processor and spans all
        nodes that hold this condition.<br>
        Also, setting a condition to the discretization set the Filled() flag
        to false.
        
  \note
  
  \param name : Name of condition
  \param cond : The Condition class
  
  \note this takes ownership of cond
  
  \warning If a condition with the exact same name already exists, it will
           NOT be overwritten but stored twice in the element
  
  */
  virtual void SetCondition(const string& name,RefCountPtr<Condition> cond);

  /*!
  \brief Get all conditions with a certain name
  
  Get all conditions with a certain name. A vector of ptrs to all conditions
  with name name is returned in out. The number of conditions found with name
  name is out.size(). out.size() is 0 if no condition with that name is found.
  
  \note Conditions attached to the discretization have to be 
        completely redundant meaning that nodal cloud in the
        condition is the same on each processor and spans all
        nodes that hold this condition

  \param name (in): Name of condition
  \param out  (out): vector of pointers to all conditions with that name
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual void GetCondition(const string& name, vector<DRT::Condition*>& out);

  /*!
  \brief Get a condition with a certain name
  
  Returns the first condition with name name found in the multimap. 
  If multiple conditions with the same name exist, the first condition is
  returned and behaviour is therefore non-deterministic. This method should
  therefore only be used in cases where the user is sure that name is unique.
  
  \note Conditions attached to the discretization have to be 
        completely redundant meaning that nodal cloud in the
        condition is the same on each processor and spans all
        nodes that hold this condition

  \param name (in): Name of condition
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual DRT::Condition* GetCondition(const string& name);

  /*!
  \brief Redistribute the discretization according to provided maps
         (Filled()==true prerequisite)
  
  Steps taken in this method are:<br>
  - build element maps (row and column)
  - do export of row/col nodes and row/col elements
  - call Fillcomplete()
  
  \note Filled()==true is a prerequisite, Filled()==true on exit
  */
  virtual void Redistribute(const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap);

  /*!
  \brief Export the nodes to a different parallel layout
         (Filled()==true NOT prerequisite)
  
  The discretization has a parallel layout of nodes reflected in NodeRowMap().
  This method communicates the nodes such that after
  the export the nodes are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - Filled()=false on exit. This means FillComplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported node distribution still matches the element distribution. A call
    to FillComplete() might therefore not be possible until elements are
    distributed accordingly.<br>
  - newmap has to be a non-overlapping map (will be tested) as
    the ownership of an exported node changes to the receiving proc<br>
  - All ghosted nodes on all processors will be destroyed.
    Ghosting has to be recreated afterwards<br>

  \note Sets Filled()=false and deletes noderowmap_ and nodecolmap_
  */
  virtual void ExportRowNodes(const Epetra_Map& newmap);


  /*!
  \brief Export overlap of nodes
        (Filled()==true NOT prerequisite)
  
  The discretization has a parallel layout of nodes reflected in NodeRowMap().
  This method communicates the nodes such that after
  the export the nodes are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - All existing ghosted nodes on all processors will be destroyed
    before the communication.<br>
  - newmap must contain all nodes of noderowmap_ (will be tested) because
    otherwise a node is shipped to a different proc and deleted from the
    originating proc. It then merely exists as a ghost node on the receiving proc.<br>
  - Filled()=false on exit. This means FillComplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported node distribution still matches the element distribution. A call
    to FillComplete() might therefore not be possible until elements are
    distributed accordingly.<br>
  - newmap should be an overlapping map.<br>
  - The ownership of an exported node does not change on the receiving proc.<br>
    The received node becomes a ghost node.

  \note Sets Filled()=false and deletes noderowmap_ and nodecolmap_
  */
  virtual void ExportColumnNodes(const Epetra_Map& newmap);

  /*!
  \brief Export the elements to a different parallel row layout
        (Filled()==true NOT prerequisite)
  
  The discretization has a parallel layout of elements reflected in ElementRowMap().
  This method communicates the elements in this discretization such that after
  the export the elements are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - Filled()=false on exit. This means FillComplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported element distribution still matches the node distribution. A call
    to FillComplete() might therefore not be possible until nodes are
    distributed accordingly.<br>
  - newmap has to be a non-overlapping map (will be tested) as
    the ownership of an exported element changes to the receiving proc<br>
  - All ghosted elements on all processors will be destroyed.
    Ghosting has to be recreated afterwards<br>

  \note Sets Filled()=false and deletes elerowmap_ and elecolmap_
  */
  virtual void ExportRowElements(const Epetra_Map& newmap);

  /*!
  \brief Export overlap of elements
        (Filled()==true NOT prerequisite)
  
  The discretization has a parallel layout of elements reflected in ElementRowMap().
  This method communicates the elements such that after
  the export the elements are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - All existing ghosted elements on all processors will be destroyed
    before the communication.<br>
  - newmap must contain all elements of elerowmap_ (will be tested) because
    otherwise an element is shipped to a different proc and deleted from the
    originating proc. It then merely exists as a ghost element on the receiving proc.<br>
  - Filled()=false on exit. This means FillComplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported element distribution still matches the node distribution. A call
    to FillComplete() might therefore not be possible until elements are
    distributed accordingly.<br>
  - newmap should be an overlapping map<br>
  - The ownership of an exported element does not change on the receiving proc.<br>
    The received element becomes a ghost element.

  \note Sets Filled()=false and deletes elerowmap_ and elecolmap_
  */
  virtual void ExportColumnElements(const Epetra_Map& newmap);

  /*!
  \brief Build nodal graph of discretization (Filled()==true prerequisite)
  
  Build a nodal graph of the discretization in parallel.<br>
  The graph has a row map of NodeRowMap().<br>
  The graph is build from elements stored on each proc, nodes are
  not referenced.<br>
  If a proc stores the appropiate ghosted elements the resulting graph
  will be the complete graph of the distributed discretization.<br>
  If procs do not store appropiate ghosted elements, the resulting
  graph is decoupled or partially decoupled among procs.<br>
  This might also lead to an unsymmetric graph. 
  
  \note Filled()=true is a prerequisite
  
  \return Graph of discretization distributed across processors according to
          the discretization distribution
  */
  virtual RefCountPtr<Epetra_CrsGraph> BuildNodeGraph() const;


  /*!
  \brief Set Ids of design nodes in nodes (Filled()==true NOT prerequisite)
  
  On exit, each node in this that is positioned on a design entity 
  will store the global design node number that it concides to.
  It will also store a flag indicating on what type of design entitty its on.
  
  \param type   : type of design entity (on_none, on_dnode, on_dline, on_dsurface, on_dvolume) that is set
  \param nfenode: vector of length ndnode, nfenode[i] is number of fe nodes on design entity i
  \param fenode: vector of length ndnode, fenode[i] is vector of global fe node ids that are on design entity i

  \note this is an individual call
  */
  virtual void SetDesignEntityIds(Node::OnDesignEntity type, const vector<int>& nfenode, const vector<vector<int> >& fenode);


  //@}

  // @{ \name Evaluate methods

  /*!
  \brief Set a reference to a data vector
  
  Using this method, a reference to a state vector can 
  be supplied to the discretization. The elements can access
  this vector by using the name of that vector.
  If the vector is supplied in DofColMap() a reference to it will be stored.
  If the vector is NOT supplied in DofColMap() a vector with column map
  is allocated and the supplied vector is exported to it.
  Everything is stored/referenced using RefCountPtr.
  
  \param name (in): Name of solution state
  \param state (in): vector of some data
  
  \note This class will not take ownership or in any way modify the solution vector.
  */
  virtual void SetState(const string& name,RefCountPtr<const Epetra_Vector> state);

  /*!
  \brief Get a reference to a data vector
  
  Providing a name of a solution state, get a reference to the solution vector
  
  \param name (in): Name of solution state
  
  \return Reference to solution state or Teuchos::null if not found
  */
  virtual RefCountPtr<const Epetra_Vector> GetState(const string& name) const
  {
    map<string,RefCountPtr<const Epetra_Vector> >::const_iterator fool = state_.find(name);
    if (fool==state_.end()) dserror("Cannot find state %s",&name[0]);
    else return fool->second;
    return null;
  }

  /*!
  \brief Clear solution state references
    
  The method deletes all references to any solution data
  */
  virtual void ClearState()
  { state_.clear(); return; }

  /*!
  \brief Call elements to evaluate
    
  */
  virtual void Evaluate(ParameterList&                params, 
                        RefCountPtr<Epetra_CrsMatrix> systemmatrix1,
                        RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                        RefCountPtr<Epetra_Vector>    systemvector1,
                        RefCountPtr<Epetra_Vector>    systemvector2,
                        RefCountPtr<Epetra_Vector>    systemvector3);


  /*!
  \brief Evaluate Neumann boundary conditions
    
  */
  virtual void EvaluateNeumann(ParameterList& params, Epetra_Vector& systemvector);

  //@}

private:

  // don't want = operator and copy constructor
  Discretization operator = (const Discretization& old);
  Discretization(const DRT::Discretization& old);

  /*!
  \brief Reset all maps and set Filled()=false (Filled()==true NOT prerequisite)
  
  \note This is a collective call
  */
  virtual void Reset();

  /*!
  \brief Build noderowmap_ (Filled()==true NOT prerequisite)
  
  Build the parallel layout of nodes in this
  discretization and store it as an Epetra_Map in noderowmap_
  noderowmap_ is non-ambigous / non-overlapping.
  It considers nodes owned by a proc only
  
  \note This is a collective call
  */
  virtual void BuildNodeRowMap();

  /*!
  \brief Build nodecolmap_ (Filled()==true NOT prerequisite)
  
  Build the parallel layout of nodes in this
  discretization and store it as an Epetra_Map in nodecolmap_
  nodecolmap_ is ambigous / potentially overlapping.
  It considers nodes owned by a proc and its ghosted nodes
  
  \note This is a collective call
  */
  virtual void BuildNodeColMap();

  /*!
  \brief Build elerowmap_ (Filled()==true NOT prerequisite)
  
  Build the parallel layout of elements in this
  discretization and store it as an Epetra_Map in elerowmap_
  elerowmap_ is non-ambigous / non-overlapping.
  It considers elements owned by a proc only
  
  \note This is a collective call
  
  */
  virtual void BuildElementRowMap();

  /*!
  \brief Build elecolmap_ (Filled()==true NOT prerequisite)
  
  Build the potentially overlapping parallel layout of elements in this
  discretization and store it as an Epetra_Map in elecolmap_
  elecolmap_ includes ghosted elements
  
  \note This is a collective call
  
  */
  virtual void BuildElementColMap();

  /*!
  \brief Build pointers elements -> Nodes (Filled()==true NOT prerequisite)
  */
  virtual void BuildElementToNodePointers();

  /*!
  \brief Build pointers Node -> Element (Filled()==true NOT prerequisite)
  
  \note This is a collective call
  */
  virtual void BuildNodeToElementPointers();

  /*!
  \brief (Re)create a map of ElementRegister objects
  
  */
  virtual void BuildElementRegister();

  /*!
  \brief Build element row and column map from nodal row and column maps
         (Filled()==true prerequisite)
  
  Create a unique element map elerowmap assigning each element an owner.
  Create an overlapping element map elecolmap representing a one layer
  overlap of element.
  Maps are created such that they match the nodal row and column maps provided.
  
  \param noderowmap (in): unique nodal row map of some distribution
  \param nodecolmap (in): overlapping nodal column map
  \param elerowmap (out): unique element row map
  \param elecolmap (out): overlapping element column map
  
  \note The provided noderowmap and nodecolmap do not need to match the
  distribution of nodes in this discretization class. Also, the output
  element maps do not match the distribution of elements in this class.
  Total numbers of nodes and elements have to match nodes and elements in 
  this class.
  
  \return Filled()==true on exit
  */
  virtual void BuildElementRowColumn(const Epetra_Map& noderowmap,
                                     const Epetra_Map& nodecolmap,
                                     RefCountPtr<Epetra_Map>& elerowmap,
                                     RefCountPtr<Epetra_Map>& elecolmap) const;

  /*!
  \brief Assign degrees of freedom to discretization (Filled()==true prerequisite)
  
  Assign nodes and elements their no. of degrees of freedom as aquired
  by Element::NumDofPerNode and Element::NumDofPerElement.
  Number degrees of freedom (dofs) ascending according to global node numbers
  followed by global element numbers
  
  \note This should become public to allow users to assign there own dofs
        in a more-than-one-discretization case
  
  \param start (int): first dof number to assign
  
  \return last dof assigned + 1
  */
  virtual int AssignDegreesOfFreedom(const int start);

  /*!
  \brief Initialize element routines
  
  Loops through the ElementRegister map elementregister_ and call the
  initialize method of each of the ElementRegister classes present.
  
  
  \note InitializeElements might be called more then once!
  
  */
  virtual void InitializeElements();

  /*!
  \brief Build the geometry for boundary conditions
  
  */
  virtual void BoundaryConditionsGeometry();

  /*!
  \brief Build the geometry of lines for a certain line condition
  
  */
  virtual void BuildLinesinCondition(const string name,RefCountPtr<DRT::Condition> cond);

  /*!
  \brief Build the geometry of surfaces for a certain surface condition
  
  */
  virtual void BuildSurfacesinCondition(const string name,RefCountPtr<DRT::Condition> cond);

  /*!
  \brief Build the geometry of volumes for a certain volume condition
  
  */
  virtual void BuildVolumesinCondition(const string name,RefCountPtr<DRT::Condition> cond);

protected:

  // internal data
  string                                        name_;            // name of this discretization
  RefCountPtr<Epetra_Comm>                      comm_;            // an Epetra_comm
  bool                                          filled_;          // flag indicating whether getOrganized has been called
  bool                                          havedof_;         // flag indicating whether degrees of freedom where assigned

  RefCountPtr<Epetra_Map>                       dofrowmap_;       // unique row map of degrees of freedom (node and element dofs))
  RefCountPtr<Epetra_Map>                       dofcolmap_;       // unique column map of degrees of freedom (node and element dofs))

  RefCountPtr<Epetra_Map>                       elerowmap_;       // unique distribution of element ownerships
  RefCountPtr<Epetra_Map>                       elecolmap_;       // distribution of elements including ghost elements
  vector<DRT::Element*>                         elerowptr_;       // vector of pointers to row elements for faster access
  vector<DRT::Element*>                         elecolptr_;       // vector of pointers to column elements for faster access
  map<int,RefCountPtr<DRT::Element> >           element_;         // map of elements
  
  RefCountPtr<Epetra_Map>                       noderowmap_;      // unique distribution of nodal ownerships
  RefCountPtr<Epetra_Map>                       nodecolmap_;      // distribution of nodes including ghost nodes
  vector<DRT::Node*>                            noderowptr_;      // vector of pointers to row nodes for faster access
  vector<DRT::Node*>                            nodecolptr_;      // vector of pointers to column nodes for faster access
  map<int,RefCountPtr<DRT::Node> >              node_;            // map of nodes

  map<string,RefCountPtr<const Epetra_Vector> > state_;           // map of references to solution states
  map<int,RefCountPtr<ElementRegister> >        elementregister_; // register for types of elements in this discretization

  multimap<string,RefCountPtr<Condition> >      condition_;       // some conditions e.g. BCs

}; // class Discretization
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Discretization& dis);









#endif  // #ifndef DISCRET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
