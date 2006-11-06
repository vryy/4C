/*!----------------------------------------------------------------------
\file node.H

\class CCADISCRETIZATION::Node

\brief A virtual class for a node

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef NODE_H
#define NODE_H


#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

#include "parobject.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

// forward declarations
class Element;

/*!
\class Node

\brief A virtual class all nodes that are used in CCADISCRETIZATION have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Node : public CCADISCRETIZATION::ParObject
{
public:

  /*!
  \brief The Discretization is a friend of Element
  */
  friend class Discretization;

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id    : A globally unique node id
  \param coords: vector of nodal coordinates, length 3
  */
  Node(int id, const double* coords);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node
  
  */
  Node(const CCADISCRETIZATION::Node& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual CCADISCRETIZATION::Node* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~Node();

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this node
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual bool Unpack(const char* data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Return global id
  */
  inline int Id() const { return id_; }


  /*!
  \brief Return coordinates (length 3)
  */
  inline const double* X() const { return x_; }

  /*!
  \brief Print this element
  */
  virtual void Print() const ;

  //@}

  // @{ \name Construction


  //@}

private:

  /*!
  \brief Set a pointer to an element
  */
  virtual void AddElementPtr(CCADISCRETIZATION::Element* ele) const ;

  int                                 id_;      // a unique global id
  double                              x_[3];    // nodal coords
  vector<CCADISCRETIZATION::Element*> element_; // pointers to adjacent elements

}; // class Node
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Node& node);









#endif  // #ifndef NODE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
