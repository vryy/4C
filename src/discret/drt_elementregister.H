/*!----------------------------------------------------------------------
\file drt_elementregister.H

\class DRT::ElementRegister

\brief A pure virtual element register class

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ELEMENTTYPE_H
#define ELEMENTTYPE_H


#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;
#include "drt_parobject.H"
#include "drt_element.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations

/*!
\class ElementRegister

\brief A pure virtual element type class

\author gee (gee@lnm.mw.tum.de)
*/
class ElementRegister : public DRT::ParObject
{
public:

  
  /*!
  \brief The Discretization is a friend of Element
  */
//  friend class Discretization;
  

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  ElementRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy
  
  */
  ElementRegister(const DRT::ElementRegister& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual DRT::ElementRegister* Clone() const = 0;

  /*!
  \brief Destructor

  */
  virtual ~ElementRegister();

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  //virtual const char* Pack(int& size) const;
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  //virtual bool Unpack(const char* data);
  virtual void Unpack(const vector<char>& data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Get type of element
  */
  inline DRT::Element::ElementType Type() const { return etype_; }

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const;

  //@}

  // @{ \name Construction

  /*!
  \brief Call Element init routines
  
  The method performs all steps needed to init the element routines
  */
  virtual int Initialize(DRT::Discretization& dis) = 0;

  //@}



private:

  DRT::Element::ElementType           etype_;     // type of element

}; // class ElementRegister
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::ElementRegister& ele);









#endif  // #ifndef ELEMENTTYPE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
