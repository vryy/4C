/*!----------------------------------------------------------------------
\file drt_node.H

\brief Base node class

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef NODE_H
#define NODE_H


#include "Teuchos_RCP.hpp"

#include "drt_parobject.H"
#include "drt_condition.H"
#include "drt_parobjectfactory.H"


namespace DRT
{

// forward declarations
class Element;
class Discretization;
namespace MESHFREE
{
  class MeshfreeDiscretization;
}

class NodeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "NodeType"; }

  static NodeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

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
  Node(int id, const double* coords, const int owner, const bool iscosserat=false);

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
  virtual void Unpack(const std::vector<char>& data);

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
  \brief Return coordinates (length 3 or 6 for a Cosserat node)
  */
  virtual inline const double* X() const { return &x_[0]; }

  /*!
  \brief return spatial dimension of node coordinates
  */
  virtual inline int Dim() const { return x_.size(); }

  /*!
  \brief Return bool whether node is a cosserat node and therefore has 6 coordinates
  */
  virtual inline bool IsCosserat() const { return (x_.size()==(unsigned)6); }

  /*!
  \brief Return processor-local number of elements adjacent to this node
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
  virtual void Print(std::ostream& os) const ;

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
        Teuchos::RCP, one can not say who actually owns the underlying object.
        The node does not communicate any conditions through Pack/Unpack,
        Conditions are therefore more of a reference here that will be
        recreated after communications of nodes have been done.

  \warning If a condition with the exact same name already exists, it will
           NOT be overwritten but stored twice in the element

  */
  virtual void SetCondition(const std::string& name,Teuchos::RCP<Condition> cond)
  { condition_.insert(std::pair<std::string,Teuchos::RCP<Condition> >(name,cond)); return; }

  /*!
  \brief Get all conditions with a certain name

  Get all conditions with a certain name. A vector of ptrs to all conditions
  with name name is returned in out. The number of conditions found with name
  name is out.size(). out.size() is 0 if no condition with that name is found.

  \param name (in): Name of condition
  \param out  (out): vector of pointers to all conditions with that name

  */
  virtual void GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const;

  /*!
  \brief Get a condition with a certain name

  Returns the first condition with name name found in the multimap.
  If multiple conditions with the same name exist, the first condition is
  returned and behaviour is therefore non-deterministic. This method should
  therefore only be used in cases where the user is sure that name is unique.

  \param name (in): Name of condition

  \return Returns NULL if condition with that name does not exist
  */
  virtual DRT::Condition* GetCondition(const std::string& name) const;

  /*!
  \brief Delete all conditions set to this node
  */
  virtual void ClearConditions() { condition_.clear(); return; }

  /*!
  \brief Change reference position by adding input vector to position
  */
  virtual void ChangePos(std::vector<double> nvector);

  /*!
  \brief Change reference position by setting input vector to position
  */
  virtual void SetPos(std::vector<double> nvector);

  //@}

  /*! \brief Query names of node data to be visualized using BINIO
   *
   *  This method is to be overloaded by a derived class.
   *  The node is supposed to fill the provided map with key names of
   *  visualization data the node wants to visualize.
   *
   *  \return On return, the derived class has filled names with key names of
   *  data it wants to visualize and with int dimensions of that data.
   */
  virtual void VisNames(std::map<std::string,int>& names) { return; }

  /*! \brief Visualize the owner of the node using BINIO
   *
   *  \param names (out): Owner is added to the key names
   */
  virtual void VisOwner(std::map<std::string,int>& names)
  {
    names.insert(std::pair<std::string,int>("Nodeowner",1));
    //names.insert(std::pair<string,int>("EleGId",1));
    return; }

  /*! \brief Query data to be visualized using BINIO of a given name
   *
   *  This method is to be overloaded by a derived class.
   *  The derived method is supposed to call this base method to visualize the
   *  owner of the node.
   *  If the derived method recognizes a supported data name, it shall fill it
   *  with corresponding data.
   *  If it does NOT recognizes the name, it shall do nothing.
   *
   *  \warning The method must not change size of variable data
   *
   *  \param name (in): Name of data that is currently processed for visualization
   *  \param data (out): data to be filled by element if it recognizes the name
   */
  virtual bool VisData(const std::string& name, std::vector<double>& data);

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
  std::vector<double>                      x_;
  //! pointers to adjacent elements
  std::vector<Element*>                    element_;
  //! some conditions e.g. BCs
  std::multimap<std::string,Teuchos::RCP<Condition> > condition_;

}; // class Node
} // namespace DRT


// << operator
std::ostream& operator << (std::ostream& os, const DRT::Node& node);









#endif  // #ifndef NODE_H
