// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_DOFSET_TRANSPARENT_INDEPENDENT_HPP
#define FOUR_C_FEM_DOFSET_TRANSPARENT_INDEPENDENT_HPP


#include "4C_config.hpp"

#include "4C_fem_dofset_independent.hpp"
#include "4C_fem_dofset_transparent.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE


namespace Cut
{
  class CutWizard;
}

namespace Core::DOFSets
{
  /// Alias dofset that shares dof numbers with another dofset
  /*!
  A special set of degrees of freedom, implemented in order to assign the same degrees of freedom to
  nodes belonging to two discretizations. This way two discretizations can assemble into the same
  position of the system matrix. As internal variable it holds a source discretization
  (Constructor). If such a nodeset is assigned to a sub-discretization, its dofs are assigned
  according to the dofs of the source.

  */
  class TransparentIndependentDofSet : public IndependentDofSet, public TransparentDofSet
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit TransparentIndependentDofSet(
        std::shared_ptr<Core::FE::Discretization> sourcedis, bool parallel);



    /// create a copy of this object
    std::shared_ptr<DofSet> clone() override { return std::make_shared<IndependentDofSet>(*this); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    int assign_degrees_of_freedom(
        const Core::FE::Discretization& dis, const unsigned dspos, const int start) override;

   protected:
    int num_dof_per_node(const Core::Nodes::Node& node) const override;


  };  // class TransparentDofSet
}  // namespace Core::DOFSets

FOUR_C_NAMESPACE_CLOSE

#endif
