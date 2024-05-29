/*---------------------------------------------------------------------*/
/*! \file

\brief Base node class

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_NODE_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_NODE_HPP


#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_discretization_condition.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::Elements
{
  class Element;
}

namespace DRT
{
  // forward declarations
  class Discretization;
}  // namespace DRT

namespace CORE::Nodes
{
  class NodeType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "NodeType"; }

    static NodeType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static NodeType instance_;
  };

  /*!
  \brief A virtual class all nodes that are used in DRT have to implement

  */
  class Node : public CORE::COMM::ParObject
  {
   public:
    //! @name Enums and Friends

    /*!
    \brief The discretization is a friend of Node
    */
    friend class Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id     (in): A globally unique node id
    \param coords (in): vector of nodal coordinates
    \param owner  (in): Owner of this node.
    */
    Node(int id, const std::vector<double>& coords, int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Node

    */
    Node(const Node& old);

    Node& operator=(const Node&) = default;

    /*!
    \brief Deep copy the derived class and return pointer to it

    */
    virtual Node* Clone() const;


    /*!
    \brief Return unique ParObject id

    every class imploementing ParObject needs a unique id defined at the
    top of this file.
    */
    int UniqueParObjectId() const override { return NodeType::Instance().UniqueParObjectId(); }

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this node

    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this node

    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    /*!
    \brief Return global id
    */
    inline int Id() const { return id_; }

    /*!
    \brief Return processor local col map id
    */
    inline int LID() const { return lid_; }

    /*!
    \brief Return owner of this node
    */
    inline int Owner() const { return owner_; }

    /*!
    \brief Return coordinates vector
    */
    inline const std::vector<double>& X() const { return x_; }

    /*!
    \brief return spatial dimension of node coordinates
    */
    inline int Dim() const { return x_.size(); }

    /*!
    \brief Return processor-local number of elements adjacent to this node
    */
    inline int NumElement() const { return element_.size(); }

    /*!
    \brief Return ptr to vector of element ptrs
    */
    inline CORE::Elements::Element** Elements()
    {
      if (NumElement())
        return element_.data();
      else
        return nullptr;
    }

    /*!
    \brief Return const ptr to vector of const element ptrs
    */
    inline const CORE::Elements::Element* const* Elements() const
    {
      if (NumElement())
        return (const CORE::Elements::Element* const*)(element_.data());
      else
        return nullptr;
    }


    /*!
    \brief Print this node
    */
    virtual void Print(std::ostream& os) const;

    //@}

    //! @name Construction

    /*!
      \brief Set processor local col id
      \param lid: processor local col id
     */
    inline void SetLID(int lid) { lid_ = lid; }

    /*!
    \brief Set ownership

    \param owner: Proc owning this node
    */
    inline void SetOwner(const int owner) { owner_ = owner; }

    /*!
    \brief Set a condition with a certain name

    Store a condition with a certain name in the node. The name need not
    be unique, meaning multiple conditions with the same name can be stored.
    Conditions can then be accessed with the GetCondition methods.

    \param name : Name of condition
    \param cond : The Condition class

    \note Normally, This method would be called by the discretization to
          set references to a Condition in the nodes. As the Condition is
          Teuchos::RCP, one can not say who actually owns the underlying object.
          The node does not communicate any conditions through Pack/Unpack,
          Conditions are therefore more of a reference here that will be
          recreated after communications of nodes have been done.

    \warning If a condition with the exact same name already exists, it will
             NOT be overwritten but stored twice in the element

    */
    void SetCondition(const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond)
    {
      condition_.insert(
          std::pair<std::string, Teuchos::RCP<CORE::Conditions::Condition>>(name, cond));
    }

    /*!
    \brief Get all conditions with a certain name

    Get all conditions with a certain name. A vector of ptrs to all conditions
    with name name is returned in out. The number of conditions found with name
    name is out.size(). out.size() is 0 if no condition with that name is found.

    \param name (in): Name of condition
    \param out  (out): vector of pointers to all conditions with that name

    */
    void GetCondition(
        const std::string& name, std::vector<CORE::Conditions::Condition*>& out) const;

    /*!
    \brief Get a condition with a certain name

    Returns the first condition with name name found in the multimap.
    If multiple conditions with the same name exist, the first condition is
    returned and behaviour is therefore non-deterministic. This method should
    therefore only be used in cases where the user is sure that name is unique.

    \param name (in): Name of condition

    \return Returns nullptr if condition with that name does not exist
    */
    CORE::Conditions::Condition* GetCondition(const std::string& name) const;

    /*!
    \brief Delete all conditions set to this node
    */
    void ClearConditions() { condition_.clear(); }

    /*!
    \brief Change reference position by adding input vector to position
    */
    void ChangePos(std::vector<double> nvector);

    /*!
    \brief Change reference position by setting input vector to position
    */
    void SetPos(std::vector<double> nvector);

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
    virtual void VisNames(std::map<std::string, int>& names) { return; }

    /*! \brief Visualize the owner of the node using BINIO
     *
     *  \param names (out): Owner is added to the key names
     */
    void VisOwner(std::map<std::string, int>& names)
    {
      names.insert(std::pair<std::string, int>("Nodeowner", 1));
    }

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

    /*!
    \brief Clear vector of pointers to my elements

    */
    inline void clear_my_element_topology() { element_.clear(); }

    /*!
    \brief Add an element to my vector of pointers to elements

    Resizes the element ptr vector and adds ptr at the end of vector
    */
    inline void add_element_ptr(CORE::Elements::Element* eleptr)
    {
      const int size = element_.size();
      element_.resize(size + 1);
      element_[size] = eleptr;
    }

   protected:
    //! a unique global id
    int id_;
    //! local col map id
    int lid_;
    //! proc owning this node
    int owner_;
    //! nodal coords
    std::vector<double> x_;
    //! pointers to adjacent elements
    std::vector<CORE::Elements::Element*> element_;
    //! some conditions e.g. BCs
    std::multimap<std::string, Teuchos::RCP<CORE::Conditions::Condition>> condition_;

  };  // class Node
}  // namespace CORE::Nodes


// << operator
std::ostream& operator<<(std::ostream& os, const CORE::Nodes::Node& node);



FOUR_C_NAMESPACE_CLOSE

#endif
