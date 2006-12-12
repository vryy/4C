/*!----------------------------------------------------------------------
\file node.H

\class DRT::Node

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

#include "drt_parobject.H"
#include "drt_dofset.H"
#include "drt_condition.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Element;

/*!
\class Node

\brief A virtual class all nodes that are used in DRT have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Node : public DRT::ParObject
{
public:

  /*!
  \brief The Discretization is a friend of Node
  */
  friend class Discretization;

  // @{ \name Constructors and destructors

  /*!
  \brief Type of design entity this node is on
         
  */
  enum OnDesignEntity 
  {
    on_none,
    on_dnode,
    on_dline,
    on_dsurface,
    on_dvolume,
  };

  /*!
  \brief Standard Constructor
  
  \param id    : A globally unique node id
  \param coords: vector of nodal coordinates, length 3
  */
  Node(int id, const double* coords, const int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node
  
  */
  Node(const DRT::Node& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual DRT::Node* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~Node();

  /*!
  \brief Return unique ParObject id
  
  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Node;}

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
  virtual inline int Id() const { return id_; }

  /*!
  \brief Return owner of this node
  */
  virtual inline int Owner() const { return owner_; }

  /*!
  \brief Return coordinates (length 3)
  */
  virtual inline const double* X() const { return x_; }

  /*!
  \brief Return number of elements adjacent to this node
  */
  virtual inline int NumElement() const { return element_.size(); }

  /*!
  \brief Return ptr to vector of element ptrs
  */
  virtual DRT::Element** Elements() 
  { if (NumElement()) return &element_[0]; else return NULL; }

  /*!
  \brief Return reference to the set of degrees of freedom
  */
  virtual inline DRT::DofSet& Dof() { return dofset_; }

  /*!
  \brief Return reference to the set of degrees of freedom
  */
  virtual inline const DRT::DofSet& Dof() const { return dofset_; }

  /*!
  \brief Get a condition with a certain name
  
  \param name : Name of condition
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual Condition* GetCondition(const string& name);

  /*!
  \brief Get the type of design entity this node is on

  */
  virtual inline DRT::Node::OnDesignEntity GetDesignEntityType() {return dentitytype_;}

  /*!
  \brief Get the id of design entity this node is on

  */
  virtual inline int GetDesignEntityId() {return dentityid_;}

  /*!
  \brief Print this node
  */
  virtual void Print(ostream& os) const ;

  //@}

  // @{ \name Construction


  /*!
  \brief Set ownership
  
  \param owner: Proc owning this node
  */
  virtual inline void SetOwner(const int owner) { owner_ = owner; return; }

  /*!
  \brief Set a design relationship
  
  \param type: Type of design element this node is located on
  \param gid: Global id of the design object this node is located on
  */
  virtual void SetDesignEntity(const OnDesignEntity type, const int gid)
  { dentitytype_ = type; dentityid_ = gid; return; }

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
  \brief Clear vector of pointers to my elements
  
  */
  virtual inline void ClearMyElementTopology() { element_.clear(); return;}

  /*!
  \brief Add an element to my vector of pointers to elements
  
  Resizes the element ptr vector and adds ptr at the end of vector
  */
  virtual void AddElementPtr(DRT::Element* eleptr) 
  { const int size = element_.size(); element_.resize(size+1); 
    element_[size] = eleptr; return;}

protected:

  int                                 id_;          // a unique global id
  int                                 owner_;       // proc owning this node
  double                              x_[3];        // nodal coords
  DRT::DofSet                         dofset_;      // a set of degrees of freedom
  vector<Element*>                    element_;     // pointers to adjacent elements
  OnDesignEntity                      dentitytype_; // type of design object I'm on
  int                                 dentityid_;   // global id of design object I'm on
  map<string,RefCountPtr<Condition> > condition_;   // some conditions e.g. BCs

}; // class Node
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Node& node);









#endif  // #ifndef NODE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
