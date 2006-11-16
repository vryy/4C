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
  \brief Get map associated with the distribution of the elements
  
  Returns NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* ElementMap() const { if (Filled()) return elemap_.get();
                                                 else dserror("FillComplete() must e called before call to ElementMap()"); 
                                                 return NULL; }

  /*!
  \brief Get map associated with the distribution of the nodes

  Returns NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeMap() const { if (Filled()) return nodemap_.get();
                                              else dserror("FillComplete() must e called before call to NodeMap()"); 
                                              return NULL; }

  /*!
  \brief Get flag indicating whether FillComplete() has been called
  */
  virtual inline bool Filled() const {return filled_;}

  /*!
  \brief Get processor local number of elements
  */
  virtual inline int NumMyElements() const { if (Filled()) return elemap_->NumMyElements();
                                             else          return element_.size();}
  
  /*!
  \brief Get global number of elements
  
  This is a collective call
  */
  virtual int NumGlobalElements() const;
  
  /*!
  \brief Get processor local number of nodes
  */
  virtual inline int NumMyNodes() const { if (Filled()) return nodemap_->NumMyElements();
                                          else          return node_.size();}
  
  /*!
  \brief Get global number of nodes
  
  This is a collective call
  */
  virtual int NumGlobalNodes() const;
  
  /*!
  \brief Get the element with global id gid
  
  This is an individual call
  
  \return Adress of element if element is owned by calling proc, returns NULL
          otherwise
  */
  virtual CCADISCRETIZATION::Element* gElement(int gid) const;

  /*!
  \brief Get the element with local id lid
  
  This is an individual call and Filled()=true is a prerequisite
  
  \return Adress of element if element is owned by calling proc, returns NULL
          otherwise
  */
  virtual CCADISCRETIZATION::Element* lElement(int lid) const;

  /*!
  \brief Get the node with global id gid
  
  This is an individual call
  
  \return Adress of node if node is owned by calling proc, returns NULL
          otherwise
  */
  virtual CCADISCRETIZATION::Node* gNode(int gid) const;

  /*!
  \brief Get the node with local id lid
  
  This is an individual call, Filled()=true is a prerequisite
  
  \return Adress of node if node is owned by calling proc, returns NULL
          otherwise
  */
  virtual CCADISCRETIZATION::Node* lNode(int lid) const;

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
  
  The discretization has a parallel layout of nodes reflected in NodeMap().
  This method communicates the nodes in this discretization such that after
  the export the nodes are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - Filled()=false on exit. This means Fillcomplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported node distribution still matches the element distribution. A call
    to FillComplete() might therefore not be possible until elements are
    distributed accordingly.
    
  \note This is a collective call
  */
  virtual void ExportNodes(const Epetra_Map& newmap);

  /*!
  \brief Export the elements to a different parallel layout
  
  The discretization has a parallel layout of elements reflected in ElementMap().
  This method communicates the elements in this discretization such that after
  the export the elements are stored as provided in newmap.<br>
  There are some important aspects to this method:<br>
  - Filled()=false on exit. This means Fillcomplete() needs to be called again.<br>
  - This is a dull export meaning that there is no notion of whether the
    exported element distribution still matches the node distribution. A call
    to FillComplete() might therefore not be possible until nodes are
    distributed accordingly.
    
  \note This is a collective call
  */
  virtual void ExportElements(const Epetra_Map& newmap);

  /*!
  \brief Set Ids of design nodes in nodes
  
  On exit, each node in this that is positioned on a design node 
  will store the global design node number that it concides to.
  It will also store a flag indicating on what type of design entitty its on.
  
  \param type   : type of design entity (on_none, on_dnode, on_dline, on_dsurface, on_dvolume) that is set
  \param nfenode: vector of length ndnode, nfenode[i] is number of fe nodes on design entity i
  \param fenode: vector of length ndnode, fenode[i] is vector of global fe node ids that are on design entity i

  \note this is an individual call
  */
  virtual void SetDesignEntityIds(Node::OnDesignEntity type, const vector<int>& nfenode, const vector<vector<int> >& fenode);

  //@}

protected:

  // don't want = operator
  Discretization operator = (const Discretization& old);

  /*!
  \brief Build nodemap_
  
  Build the potentially overlapping parallel layout of nodes in this
  discretization and store it as an Epetra_Map in nodemap_
  
  \note This is a collective call
  
  */
  virtual void BuildNodeMap();

  /*!
  \brief Build elemap_
  
  Build the potentially overlapping parallel layout of elements in this
  discretization and store it as an Epetra_Map in elemap_
  
  \note This is a collective call
  
  */
  virtual void BuildElementMap();

  /*!
  \brief Build pointers elements -> Nodes
  */
  virtual void BuildElementToNodePointers();

  /*!
  \brief Build pointers Node -> Element
  
  \note This is a collective call
  */
  virtual void BuildNodeToElementPointers();

  // internal data
  RefCountPtr<Epetra_Comm>                          comm_;        // an Epetra_comm
  bool                                              filled_;      // flag indicating whether getOrganized has been called

  RefCountPtr<Epetra_Map>                           elemap_;      // indicates the potentially overlapping element distribution
  map<int,RefCountPtr<CCADISCRETIZATION::Element> > element_;     // map of elements
  RefCountPtr<Epetra_Map>                           nodemap_;     // indicates the potentially overlapping node distribution
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >    node_;        // map of nodes


}; // class Discretization
} // namespace CCADISCRET


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Discretization& dis);









#endif  // #ifndef DISCRET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
