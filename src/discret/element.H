/*!----------------------------------------------------------------------
\file element.H

\class CCADISCRETIZATION::Shell8

\brief A C++ wrapper for the shell8 element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ELEMENT_H
#define ELEMENT_H


#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;
#include "parobject.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Element

\brief A virtual class all elements that are used in CCADISCRETIZATION have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Element : public CCADISCRETIZATION::ParObject
{
public:

  /*!
  \brief Type of element
         
   \param element_none : default value
   \param element_shell8 :shell8 element
  */
  enum ElementType 
  {
    element_none,
    element_shell8
  };


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id   : A globally unique element id
  \param etype: Type of element
  */
  Element(int id, ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Element(const CCADISCRETIZATION::Element& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual CCADISCRETIZATION::Element* Clone() const = 0;

  /*!
  \brief Destructor

  */
  virtual ~Element();

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Return global id of this element
  */
  inline int Id() const { return id_; }

  /*!
  \brief Get type of element
  */
  inline ElementType Type() const { return etype_; }

  /*!
  \brief Return number of nodes of this element
  */
  inline int NumNode() const { return node_.size(); }

  /*!
  \brief Print this element (pure virtual)
  */
  virtual void Print() const = 0;

  //@}

  // @{ \name Construction

  /*!
  \brief Read input for this element
  */
  virtual bool ReadElement();

  /*!
  \brief Set a list of node ids this element is connected to
  
  Sets the nodal ids of the nodes adjacent to this element and the number
  of nodes
  
  \param nnode : number of nodes
  \param nodes : list of unique global nodal ids
  
  */
  void SetNodes(const int nnode, const int* nodes);

  //@}

private:

  int         id_;       // a unique global element id
  ElementType etype_;    // type of element
  vector<int> node_;     // list of my nodal ids

}; // class Element
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Element& seg);









#endif  // #ifndef ELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
