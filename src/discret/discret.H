/*!----------------------------------------------------------------------
\file discret.H

\class CCADISCRETIZATION::Discretization

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

#include "element.H"
#include "node.H"
#include "dserror.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Discretization

\brief A class to manage a discretization in parallel

The \ref CCADISCRET::Discretization class supports the ostream& operator <<

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
  Discretization(RefCountPtr<Epetra_Comm> comm);


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
  \brief Get map associated with the distribution of the ownership of elements
  
  This map includes all elements stored on this proc and also owned by this proc.
  This map is non-ambiguous, meaning that it is a non-overlapping map.
  
  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* ElementRowMap() const { if (Filled()) return elerowmap_.get();
                                                    else dserror("FillComplete() must be called before call to ElementRowMap()"); 
                                                    return NULL; }
  /*!
  \brief Get map associated with the distribution of elements including ghosted elements
  
  This map includes all elements stored on this proc including any ghosted elements
  This map is ambiguous, meaning that it is an overlapping map
  
  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* ElementColMap() const { if (Filled()) return elecolmap_.get();
                                                    else dserror("FillComplete() must be called before call to ElementColMap()"); 
                                                    return NULL; }

  /*!
  \brief Get map associated with the distribution of the ownership of nodes

  This map includes all nodes stored on this proc and also owned by this proc.
  This map is non-ambiguous, meaning that it is a non-overlapping map.

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeRowMap() const { if (Filled()) return noderowmap_.get();
                                                 else dserror("FillComplete() must be called before call to NodeRowMap()"); 
                                                 return NULL; }

  /*!
  \brief Get map associated with the distribution of nodes including ghosted nodes

  This map includes all nodes stored on this proc including any ghosted nodes
  This map is ambiguous, meaning that it is an overlapping map

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeColMap() const { if (Filled()) return nodecolmap_.get();
                                                 else dserror("FillComplete() must be called before call to NodeRowMap()"); 
                                                 return NULL; }

  /*!
  \brief Get global number of elements (true number of total elements)
  
  This is a collective call
  */
  virtual int NumGlobalElements() const { if (Filled()) return ElementRowMap()->NumGlobalElements();
                                          else dserror("FillComplete() must be called before call to NumGlobalElements()");
                                          return -1; } 
  
  /*!
  \brief Get processor local number of elements owned by this processor
  */
  virtual int NumMyRowElements() const { if (Filled()) return ElementRowMap()->NumMyElements();
                                         else dserror("FillComplete() must be called before call to NumMyRowElements()"); 
                                         return -1;}

  /*!
  \brief Get processor local number of elements including ghost elements
  */
  virtual int NumMyColElements() const { if (Filled()) return ElementColMap()->NumMyElements();
                                         else return (int)element_.size(); 
                                         return -1;}
  
  /*!
  \brief Get global number of nodes (true number of total nodes without ghosting)
  
  This is a collective call
  */
  virtual int NumGlobalNodes() const { if (Filled()) return NodeRowMap()->NumGlobalElements();
                                       else dserror("FillComplete() must be called before call to NumGlobalNodes()");
                                       return -1; }

  /*!
  \brief Get processor local number of nodes owned by this processor
  */
  virtual int NumMyRowNodes() const { if (Filled()) return NodeRowMap()->NumMyElements();
                                      else dserror("FillComplete() must be called before call to NumMyRowNodes()");
                                      return -1; }
  
  /*!
  \brief Get processor local number of nodes including ghost nodes
  */
  virtual int NumMyColNodes() const { if (Filled()) return NodeColMap()->NumMyElements();
                                      else return (int)node_.size();
                                      return -1; }
  
  
  /*!
  \brief Get the element with global id gid
  
  Returns the element with global row id gid if element is on this proc.
  Will return row or column element, ghosted or not.
  This is an individual call
  
  \return Adress of element if element is owned by calling proc, returns NULL
          otherwise
  */
  virtual CCADISCRETIZATION::Element* gElement(int gid) const;

  /*!
  \brief Get the element with local row id lid
  
  Returns the element with local row index lid.
  Will not return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of element if element is owned by calling proc
  */
  virtual CCADISCRETIZATION::Element* lRowElement(int lid) const;

  /*!
  \brief Get the element with local column id lid
  
  Returns the element with local column index lid.
  Will also return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of element if element is stored by calling proc
  */
  virtual CCADISCRETIZATION::Element* lColElement(int lid) const;

  /*!
  \brief Get the node with global row id gid
  
  Returns the node with global row id gid if node is on this proc.
  Will return row or column node, ghosted or  not.
  This is an individual call
  
  \return Adress of node if node is stored on calling proc
  */
  virtual CCADISCRETIZATION::Node* gNode(int gid) const;

  /*!
  \brief Get the node with local row id lid
  
  Returns the node with local row index lid.
  Will not return any ghosted node.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of node if node is owned and stored by calling proc
  */
  virtual CCADISCRETIZATION::Node* lRowNode(int lid) const;

  /*!
  \brief Get the node with local column id lid
  
  Returns the node with local column index lid.
  Will return any node stored on this proc.
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of node if node is stored by calling proc
  */
  virtual CCADISCRETIZATION::Node* lColNode(int lid) const;

  /*!
  \brief Print this discretization to stdout
  
  \note This is a collective call
  */
  virtual void Print(ostream& os) const;

  //@}

  // @{ \name Construction methods

  /*!
  \brief Set an epetra_Commobject which can be either 
         Epetra_SerialComm or Epetra_MpiComm
  
  \note Sets the Filled() flag to false, Fillcomplete must  be called 
        afterwards to rebuild node and element map
  */
  virtual void SetComm(RefCountPtr<Epetra_Comm> comm)
  { filled_=false; comm_ = comm; return; }

  /*!
  \brief Add an element to the discretization
  
  The discretization takes ownership of the added element.
  Note that if an element with the same Id() exists, it will be 
  deleted and replaced by the new one.
  */
  virtual void AddElement(RefCountPtr<CCADISCRETIZATION::Element> ele);

  /*!
  \brief Add a node to the discretization
  
  The discretization takes ownership of the added node.
  Note that if a node with the same Id() exists, it will be 
  deleted and replaced by the new one.
  */
  virtual void AddNode(RefCountPtr<CCADISCRETIZATION::Node> node);

  /*!
  \brief Complete construction of this discretization
  
  after adding or deleting nodes or elements or redistributing them in parallel, 
  this method has to be called to (re)construct pointer topologies
  
  \note This is a collective call
  */
  virtual int FillComplete();

  /*!
  \brief Export the nodes to a different parallel layout
  
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
  */
  virtual void ExportNodes(const Epetra_Map& newmap);

  /*!
  \brief Create ghosting of nodes
  
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
  */
  virtual void ExportGhostNodes(const Epetra_Map& newmap);

  /*!
  \brief Export the elements to a different parallel layout
  
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

    
  \note This is a collective call
  */
  virtual void ExportElements(const Epetra_Map& newmap);

  /*!
  \brief Create ghosting of elements
  
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
  */
  virtual void ExportGhostElements(const Epetra_Map& newmap);

  /*!
  \brief Build nodal graph of discretization
  
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
  \brief Build element row and column map from nodal row and column map
  */
  virtual void BuildElementGhosting(const Epetra_Map& noderowmap,
                                    const Epetra_Map& nodecolmap,
                                    RefCountPtr<Epetra_Map>& elerowmap,
                                    RefCountPtr<Epetra_Map>& elecolmap) const;

  /*!
  \brief Set Ids of design nodes in nodes
  
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

private:

  // don't want = operator and copy constructor
  Discretization operator = (const Discretization& old);
  Discretization(const CCADISCRETIZATION::Discretization& old);

  /*!
  \brief Build noderowmap_
  
  Build the parallel layout of nodes in this
  discretization and store it as an Epetra_Map in noderowmap_
  noderowmap_ is non-ambigous / non-overlapping.
  It considers nodes owned by a proc only
  
  \note This is a collective call
  */
  virtual void BuildNodeRowMap();

  /*!
  \brief Build nodecolmap_
  
  Build the parallel layout of nodes in this
  discretization and store it as an Epetra_Map in nodecolmap_
  nodecolmap_ is ambigous / potentially overlapping.
  It considers nodes owned by a proc and its ghosted nodes
  
  \note This is a collective call
  */
  virtual void BuildNodeColMap();

  /*!
  \brief Build elerowmap_
  
  Build the parallel layout of elements in this
  discretization and store it as an Epetra_Map in elerowmap_
  elerowmap_ is non-ambigous / non-overlapping.
  It considers elements owned by a proc only
  
  \note This is a collective call
  
  */
  virtual void BuildElementRowMap();

  /*!
  \brief Build elecolmap_
  
  Build the potentially overlapping parallel layout of elements in this
  discretization and store it as an Epetra_Map in elecolmap_
  elecolmap_ includes ghosted elements
  
  \note This is a collective call
  
  */
  virtual void BuildElementColMap();

  /*!
  \brief Build pointers elements -> Nodes
  */
  virtual void BuildElementToNodePointers();

  /*!
  \brief Build pointers Node -> Element
  
  \note This is a collective call
  */
  virtual void BuildNodeToElementPointers();

protected:

  // internal data
  RefCountPtr<Epetra_Comm>                          comm_;        // an Epetra_comm
  bool                                              filled_;      // flag indicating whether getOrganized has been called

  RefCountPtr<Epetra_Map>                           elerowmap_;   // unique distribution of element ownerships
  RefCountPtr<Epetra_Map>                           elecolmap_;   // distribution of elements including ghost elements
  map<int,RefCountPtr<CCADISCRETIZATION::Element> > element_;     // map of elements
  
  RefCountPtr<Epetra_Map>                           noderowmap_;  // unique distribution of nodal ownerships
  RefCountPtr<Epetra_Map>                           nodecolmap_;  // distribution of nodes including ghost nodes
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >    node_;        // map of nodes


}; // class Discretization
} // namespace CCADISCRET


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Discretization& dis);









#endif  // #ifndef DISCRET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
