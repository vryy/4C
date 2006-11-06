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
  inline const Epetra_Comm& Comm() const {return *comm_;}

  /*!
  \brief Get flag indicating whether FillComplete() has been called
  */
  inline bool Filled() const {return filled_;}

  /*!
  \brief Get processor local number of elements
  */
  inline int NumMyElements() const {return element_.size();}
  
  /*!
  \brief Get global number of elements
  
  This is a collective call
  */
  int NumGlobalElements() const;
  
  /*!
  \brief Get processor local number of nodes
  */
  inline int NumMyNodes() const {return node_.size();}
  
  /*!
  \brief Get global number of nodes
  
  This is a collective call
  */
  int NumGlobalNodes() const;
  
  /*!
  \brief Get the element with global id gid
  
  This is an individual call
  
  \return Adress of element if element is owned by calling proc, returns NULL
          otherwise
  */
  const CCADISCRETIZATION::Element* Element(int gid) const;


  /*!
  \brief Get the node with global id gid
  
  This is an individual call
  
  \return Adress of node if node is owned by calling proc, returns NULL
          otherwise
  */
  const CCADISCRETIZATION::Node* Node(int gid) const;

  /*!
  \brief Print this discretization to stdout
  
  \note This is a collective call
  */
  void Print() const;

  //@}

  // @{ \name Construction methods

  /*!
  \brief Add an element to the discretization
  
  The discretization takes ownership of the added element.
  Note that if an element with the same Id() exists, it will be 
  deleted and replaced by the new one.
  */
  void AddElement(RefCountPtr<CCADISCRETIZATION::Element> ele);

  /*!
  \brief Add a node to the discretization
  
  The discretization takes ownership of the added node.
  Note that if a node with the same Id() exists, it will be 
  deleted and replaced by the new one.
  */
  void AddNode(RefCountPtr<CCADISCRETIZATION::Node> node);

  /*!
  \brief Complete construction of this discretization
  
  after adding or deleting nodes or elements or redistributing them in parallel, 
  this method has to be called to (re)construct pointer topologies
  
  \note This is a collective call
  */
  void FillComplete();

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
  void BuildNodeMap();

  /*!
  \brief Build elemap_
  
  Build the potentially overlapping parallel layout of elements in this
  discretization and store it as an Epetra_Map in elemap_
  
  \note This is a collective call
  
  */
  void BuildElementMap();

  /*!
  \brief Build pointers elements -> Nodes
  */
  void BuildElementToNodePointers();

  /*!
  \brief Build pointers Node -> Element
  
  \note This is a collective call
  */
  void BuildNodeToElementPointers();

  RefCountPtr<Epetra_Comm>                          comm_;        // an Epetra_comm
  bool                                              filled_;      // flag indiacting whether getOrganized has been called

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
