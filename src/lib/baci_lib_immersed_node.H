/*----------------------------------------------------------------------*/
/*! \file

\brief specialized node for immersed problems.

\level 2

*----------------------------------------------------------------------*/
#ifndef BACI_LIB_IMMERSED_NODE_H
#define BACI_LIB_IMMERSED_NODE_H


#include "baci_comm_parobjectfactory.H"
#include "baci_lib_node.H"

namespace DRT
{
  class Node;

  class ImmersedNodeType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ImmersedNodeType"; }

    static ImmersedNodeType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;


   private:
    static ImmersedNodeType instance_;
  };

  /*!
  \brief A virtual class all nodes that are used in DRT have to implement

  */
  class ImmersedNode : public DRT::Node
  {
   public:
    //! @name Enums and Friends

    /*!
    \brief The Discretization is a friend of Node
    */
    friend class Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id     (in): A globally unique node id
    \param coords (in): vector of nodal coordinates, length 3
    \param owner  (in): Owner of this node.
    */
    ImmersedNode(int id, const std::vector<double>& coords, const int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Node

    */
    ImmersedNode(const DRT::ImmersedNode& old);

    /*!
    \brief Deep copy the derived class and return pointer to it

    */
    DRT::ImmersedNode* Clone() const override;


    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of this file.
    */
    int UniqueParObjectId() const override
    {
      return ImmersedNodeType::Instance().UniqueParObjectId();
    }

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

    /*!
    \brief set 'true' if node is covered by an immersed discretization
    */
    void SetIsMatched(int ismatched) { ismatched_ = ismatched; }

    /*!
    \brief is node covered by an immersed discretization ?
    */
    int IsMatched() const { return ismatched_; }

    /*!
    \brief set 'true' if parent element is cut by an immersed boundary
    */
    void SetIsBoundaryImmersed(bool isbdryimmersed) { IsBoundaryImmersed_ = isbdryimmersed; }

    /*!
    \brief set 'true' if parent element is adjacent to immersed boundary and fully covered by
    immersed body.
    */
    void SetIsPseudoBoundary(bool isbdryimmersed) { IsPseudoBoundary_ = isbdryimmersed; }

    /*!
    \brief is an boundary immersed in parent element ?
    */
    int IsBoundaryImmersed() const { return IsBoundaryImmersed_; }

    /*!
    \brief is pseudo boundary node ?
    */
    int IsPseudoBoundary() const { return IsPseudoBoundary_; }

    /*!
    \brief Print this node
    */
    void Print(std::ostream& os) const override;

    /*! \brief Query names of node data to be visualized using BINIO
     *
     */
    void VisNames(std::map<std::string, int>& names) override;

    /*! \brief Query data to be visualized using BINIO of a given name
     *
     */
    bool VisData(const std::string& name, std::vector<double>& data) override;

   protected:
    //! @name immersed information
    //@{
    int ismatched_;           //!< is covered by immersed dis?
    int IsBoundaryImmersed_;  //!< is attached to an element cut by immersed boundary?
    int IsPseudoBoundary_;    //!< is part of the pseudo-boundary between physical and artificial
                              //!< domain?
    //@}

  };  // class ImmersedNode
}  // namespace DRT

// << operator
std::ostream& operator<<(std::ostream& os, const DRT::ImmersedNode& immersednode);



#endif  // BACI_LIB_IMMERSED_NODE_H
