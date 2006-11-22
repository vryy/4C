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
#include "condition.H"

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
  \brief Get a condition with a certain name
  
  \param name : Name of condition
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual const Condition* GetCondition(const string& name) const;

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

  int                                 id_;          // a unique global id
  int                                 owner_;       // proc owning this node
  double                              x_[3];        // nodal coords
  vector<Element*>                    element_;     // pointers to adjacent elements
  OnDesignEntity                      dentitytype_; // type of design object I'm on
  int                                 dentityid_;   // global id of design object I'm on
  map<string,RefCountPtr<Condition> > condition_;   // some conditions e.g. BCs

}; // class Node
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Node& node);









#endif  // #ifndef NODE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
