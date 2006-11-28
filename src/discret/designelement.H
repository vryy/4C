/*!----------------------------------------------------------------------
\file designelement.H

\class DRT::DesignElement

\brief An element used in CAD design description

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DESIGNELEMENT_H
#define DESIGNELEMENT_H


#include "element.H"
using namespace std;
using namespace Teuchos;


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class DesignElement

\brief An element used in CAD design description

\author gee (gee@lnm.mw.tum.de)
*/
class DesignElement : public DRT::Element
{
public:

  /*!
  \brief The DesignDiscretization is a friend of DesignElement
  */
  friend class DesignDiscretization;

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id : A unique global id
  */
  DesignElement(int id, enum ElementType type, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Element
  
  */
  DesignElement(const DRT::DesignElement& old);

  /*!
  \brief Deep copy this instance of DesignElement and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::DesignElement* Clone() const;

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

  /*!
  \brief Destructor

  */
  virtual ~DesignElement();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Get number of higher entities adjacent to this element
  */
  int NumHigherEntityIds() const { return hentityid_.size(); }

  /*!
  \brief Get number of lower entities adjacent to this element
  */
  int NumLowerEntityIds() const { return lentityid_.size(); }

  /*!
  \brief Get ids of lower entities adjacent to this element
  */
  inline const int* LowerEntityIds() const 
  { if (lentityid_.size()) return &lentityid_[0]; else return NULL; }

  /*!
  \brief Get ptrs to lower entities adjacent to this element
  */
  inline Element** LowerEntities()
  { if (lentity_.size()) return &lentity_[0]; else return NULL; }

  /*!
  \brief Get ids of higher entities adjacent to this element
  */
  inline const int* HigherEntityIds() const 
  { if (hentityid_.size()) return &hentityid_[0]; else return NULL; }
  /*!
  \brief Get ptrs to higher entities adjacent to this element
  */
  inline Element** HigherEntities() 
  { if (hentity_.size()) return &hentity_[0]; else return NULL; }
  //@}

  // @{ \name Construction methods

  /*!
  \brief set ids of adjacent lower entitiy elements
  
  \param nele: number of lower entities adjacent to this element
  \param ids: unique ids of lower entities
  \param orientation: orientation flag 0 for same orientation as first entity, 
                                       1 for other orientation
  
  */
  void SetLowerEntities(const int nele, const int* ids, const int* orientation);

  /*!
  \brief set ids of adjacent higher entitiy elements
  
  \param nele: number of higher entities adjacent to this element
  \param ids: unique ids of higher entities
  \param orientation: orientation flag 0 for same orientation as first entity, 
                                       1 for other orientation
  
  */
  void SetHigherEntities(const int nele, const int* ids, const int* orientation);

  /*!
  \brief Build pointers to lower entity elements
  
  \param lower: discretization holding lower entities  
  
  */
  bool BuildLowerElementPointers(const Discretization& lower);

  //@}


private:

  // don't want = operator
  DesignElement operator = (const DesignElement& old);

  // class specific data
  vector<int>      lentityid_;    // ids of lower  entity elements
  vector<int>      lorientation_; // orientation of lower entities
  vector<Element*> lentity_;      // ptrs to lower entities
  
  vector<int>      hentityid_;    // ids of higher entity elements
  vector<int>      horientation_; // orientation of higher entities
  vector<Element*> hentity_;      // ptrs to lower entities

}; // class DesignElement
} // namespace DRT











#endif  // #ifndef DESIGNELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
