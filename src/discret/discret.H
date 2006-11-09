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
                                                 else return NULL; }

  /*!
  \brief Get map associated with the distribution of the nodes

  Returns NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* NodeMap() const { if (Filled()) return nodemap_.get();
                                              else return NULL; }

  /*!
  \brief Get flag indicating whether FillComplete() has been called
  */
  virtual inline bool Filled() const {return filled_;}

  /*!
  \brief Get processor local number of elements
  */
  virtual inline int NumMyElements() const {return element_.size();}
  
  /*!
  \brief Get global number of elements
  
  This is a collective call
  */
  virtual int NumGlobalElements() const;
  
  /*!
  \brief Get processor local number of nodes
  */
  virtual inline int NumMyNodes() const {return node_.size();}
  
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
  
  This is an individual call Filled()=true is a prerequisite
  
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
