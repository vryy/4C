/*!----------------------------------------------------------------------
\file drt_node.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef NODE_H
#define NODE_H


#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

#include "drt_parobject.H"
#include "drt_condition.H"
#include "drt_parobjectfactory.H"


namespace DRT
{

// forward declarations
class Element;
class Discretisation;
namespace MESHFREE
{
  class MeshfreeDiscretization;
}

class NodeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "NodeType"; }

  static NodeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static NodeType instance_;
};

/*!
\brief A virtual class all nodes that are used in DRT have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Node : public DRT::ParObject
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of Node
  */
  friend class Discretization;
  friend class MESHFREE::MeshfreeDiscretization;

  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id     (in): A globally unique node id
  \param coords (in): vector of nodal coordinates, length 3
  \param owner  (in): Owner of this node.
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
  virtual int UniqueParObjectId() const { return NodeType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Return global id
  */
  virtual inline int Id() const { return id_; }

  /*!
  \brief Return processor local col map id
  */
  virtual inline int LID() const { return lid_; }

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
  virtual inline DRT::Element** Elements()
  { if (NumElement()) return &element_[0]; else return NULL; }

  /*!
  \brief Return const ptr to vector of const element ptrs
  */
  virtual inline const DRT::Element*const* Elements() const
  { if (NumElement()) return (const DRT::Element*const*)(&element_[0]); else return NULL; }


  /*!
  \brief Print this node
  */
  virtual void Print(ostream& os) const ;

  //@}

  //! @name Construction

  /*!
    \brief Set processor local col id
    \param lid: processor local col id
   */
  virtual inline void SetLID(int lid) { lid_ = lid; }

  /*!
  \brief Set ownership

  \param owner: Proc owning this node
  */
  virtual inline void SetOwner(const int owner) { owner_ = owner; return; }

  /*!
  \brief Set a condition with a certain name

  Store a condition with a certain name in the node. The name need not
  be unique, meaning multiple conditions with the same name can be stored.
  Conditions can then be accessed with the GetCondition methods.

  \param name : Name of condition
  \param cond : The Condition class

  \note Normally, This method would be called by the Discretization to
        set references to a Condition in the nodes. As the Condition is
        RefCountPtr, one can not say who actually owns the underlying object.
        The node does not communicate any conditions through Pack/Unpack,
        Conditions are therefore more of a reference here that will be
        recreated after communications of nodes have been done.

  \warning If a condition with the exact same name already exists, it will
           NOT be overwritten but stored twice in the element

  */
  virtual void SetCondition(const string& name,RefCountPtr<Condition> cond)
  { condition_.insert(pair<string,RefCountPtr<Condition> >(name,cond)); return; }

  /*!
  \brief Get all conditions with a certain name

  Get all conditions with a certain name. A vector of ptrs to all conditions
  with name name is returned in out. The number of conditions found with name
  name is out.size(). out.size() is 0 if no condition with that name is found.

  \param name (in): Name of condition
  \param out  (out): vector of pointers to all conditions with that name

  */
  virtual void GetCondition(const string& name, vector<DRT::Condition*>& out) const;

  /*!
  \brief Get a condition with a certain name

  Returns the first condition with name name found in the multimap.
  If multiple conditions with the same name exist, the first condition is
  returned and behaviour is therefore non-deterministic. This method should
  therefore only be used in cases where the user is sure that name is unique.

  \param name (in): Name of condition

  \return Returns NULL if condition with that name does not exist
  */
  virtual DRT::Condition* GetCondition(const string& name) const;

  /*!
  \brief Delete all conditions set to this node
  */
  virtual void ClearConditions() { condition_.clear(); return; }

  /*!
  \brief Change position reference
  */
  virtual void ChangePos(vector<double> nvector);

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
  virtual inline void AddElementPtr(DRT::Element* eleptr)
  { const int size = element_.size(); element_.resize(size+1);
    element_[size] = eleptr; return;}

protected:

  //! a unique global id
  int                                      id_;
  //! local col map id
  int                                      lid_;
  //! proc owning this node
  int                                      owner_;
  //! nodal coords
  double                                   x_[3];
  //! pointers to adjacent elements
  vector<Element*>                         element_;
  //! some conditions e.g. BCs
  multimap<string,RefCountPtr<Condition> > condition_;

}; // class Node
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Node& node);









#endif  // #ifndef NODE_H
