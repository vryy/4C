// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_IMMERSED_NODE_HPP
#define FOUR_C_FEM_GENERAL_IMMERSED_NODE_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Nodes
{
  class Node;

  class ImmersedNodeType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "ImmersedNodeType"; }

    static ImmersedNodeType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;


   private:
    static ImmersedNodeType instance_;
  };

  /*!
  \brief A virtual class all nodes used in the discretization management module have to implement

  */
  class ImmersedNode : public Node
  {
   public:
    //! @name Enums and Friends

    /*!
    \brief The discretization is a friend of Node
    */
    friend class Core::FE::Discretization;

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
    ImmersedNode(const ImmersedNode& old);

    /*!
    \brief Deep copy the derived class and return pointer to it

    */
    ImmersedNode* clone() const override;


    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of this file.
    */
    int unique_par_object_id() const override
    {
      return ImmersedNodeType::instance().unique_par_object_id();
    }

    /*!
    \brief Pack this class so it can be communicated

    \ref pack and \ref unpack are used to communicate this node

    */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref pack and \ref unpack are used to communicate this node

    */
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    /*!
    \brief set 'true' if node is covered by an immersed discretization
    */
    void set_is_matched(int ismatched) { ismatched_ = ismatched; }

    /*!
    \brief is node covered by an immersed discretization ?
    */
    int is_matched() const { return ismatched_; }

    /*!
    \brief set 'true' if parent element is cut by an immersed boundary
    */
    void set_is_boundary_immersed(bool isbdryimmersed) { IsBoundaryImmersed_ = isbdryimmersed; }

    /*!
    \brief set 'true' if parent element is adjacent to immersed boundary and fully covered by
    immersed body.
    */
    void set_is_pseudo_boundary(bool isbdryimmersed) { IsPseudoBoundary_ = isbdryimmersed; }

    /*!
    \brief is an boundary immersed in parent element ?
    */
    int is_boundary_immersed() const { return IsBoundaryImmersed_; }

    /*!
    \brief is pseudo boundary node ?
    */
    int is_pseudo_boundary() const { return IsPseudoBoundary_; }

    /*!
    \brief Print this node
    */
    void print(std::ostream& os) const override;

    /*! \brief Query names of node data to be visualized using BINIO
     *
     */
    void vis_names(std::map<std::string, int>& names) override;

    /*! \brief Query data to be visualized using BINIO of a given name
     *
     */
    bool vis_data(const std::string& name, std::vector<double>& data) override;

   protected:
    //! @name immersed information
    //@{
    int ismatched_;           //!< is covered by immersed dis?
    int IsBoundaryImmersed_;  //!< is attached to an element cut by immersed boundary?
    int IsPseudoBoundary_;    //!< is part of the pseudo-boundary between physical and artificial
                              //!< domain?
    //@}

  };  // class ImmersedNode
}  // namespace Core::Nodes

// << operator
std::ostream& operator<<(std::ostream& os, const Core::Nodes::ImmersedNode& immersednode);



FOUR_C_NAMESPACE_CLOSE

#endif
