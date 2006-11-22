/*!----------------------------------------------------------------------
\file element.H

\class CCADISCRETIZATION::Element

\brief A pure virtual element class

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
#include "condition.H"


/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

// forward declarations
class Discretization;
class Node;

/*!
\class Element

\brief A virtual class all elements that are used in CCADISCRETIZATION have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Element : public CCADISCRETIZATION::ParObject
{
public:

  
  /*!
  \brief The Discretization is a friend of Element
  */
  friend class Discretization;
  
  /*!
  \brief Type of element
         
  */
  enum ElementType 
  {
    element_none,
    element_designline,
    element_designsurface,
    element_designvolume,
    element_shell8
  };


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id   : A globally unique element id
  \param etype: Type of element
  */
  Element(int id, ElementType etype, int owner);

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
  \brief Return owner of this element
  */
  virtual inline int Owner() const { return owner_; }

  /*!
  \brief Get type of element
  */
  inline ElementType Type() const { return etype_; }

  /*!
  \brief Return number of nodes of this element
  */
  inline int NumNode() const { return nodeid_.size(); }

  /*!
  \brief Return id's of nodes adjacent to this element
  */
  inline const int* NodeIds() const { if (nodeid_.size()) return &nodeid_[0]; 
                                      else                return NULL;     }

  /*!
  \brief Get a condition with a certain name
  
  \param name : Name of condition
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual const Condition* GetCondition(const string& name) const;

  /*!
  \brief Get vector of ptrs to nodes
  
  */
  virtual Node** Nodes()
  { if (node_.size()) return &node_[0]; else return NULL; }

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const;

  //@}

  // @{ \name Construction

  /*!
  \brief Read input for this element
  */
  virtual bool ReadElement();

  /*!
  \brief Set ownership
  
  \param owner: Proc owning this node
  */
  virtual inline void SetOwner(const int owner) { owner_ = owner; return; }

  /*!
  \brief Set a list of node ids this element is connected to
  
  Sets the nodal ids of the nodes adjacent to this element and the number
  of nodes
  
  \param nnode : number of nodes
  \param nodes : list of unique global nodal ids
  
  */
  void SetNodeIds(const int nnode, const int* nodes);

  /*!
  \brief Set a condition with a certain name
  
  \param name : Name of condition
  \param cond : The Condition class
  
  \note this takes ownership of cond
  
  \warning If a condition with the exact same name already exists, it will be
           overwritten
  
  */
  virtual void SetCondition(const string& name,RefCountPtr<Condition> cond)
  { condition_[name] = cond; return; }

  //@}

protected:

  /*!
  \brief Build pointer vector from map of nodes
  
  \return Returns true on success and false if nodes from NodeIds() cannot be found
          in map
  */
  bool BuildNodalPointers(map<int,RefCountPtr<CCADISCRETIZATION::Node> >& nodes);

  int                                 id_;        // a unique global element id
  int                                 owner_;     // owner of this element
  ElementType                         etype_;     // type of element
  vector<int>                         nodeid_;    // list of my nodal ids
  vector<CCADISCRETIZATION::Node*>    node_;      // pointers to adjacent nodes
  map<string,RefCountPtr<Condition> > condition_; // some conditions e.g. BCs

}; // class Element
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Element& ele);









#endif  // #ifndef ELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
