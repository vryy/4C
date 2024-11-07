// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_DOFSET_INDEPENDENT_HPP
#define FOUR_C_FEM_DOFSET_INDEPENDENT_HPP

#include "4C_config.hpp"

#include "4C_fem_dofset.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::DOFSets
{

  /*!
  \brief A set of degrees of freedom

  \author
  */
  class IndependentDofSet : virtual public DofSet
  {
   public:
    /*!
    \brief Standard Constructor



    create a dofset that is independent of the other dofsets


    \return void

    */
    IndependentDofSet(bool ignoreminnodegid = false);

    /*!
    \brief Copy constructor

    */
    IndependentDofSet(const IndependentDofSet& old);


    /// create a copy of this object
    std::shared_ptr<DofSet> clone() override { return std::make_shared<IndependentDofSet>(*this); }

    /// Add Dof Set to list #static_dofsets_
    void add_dof_setto_list() override;

   protected:
    /// get first number to be used as Dof GID in assign_degrees_of_freedom
    int get_first_gid_number_to_be_used(const Core::FE::Discretization& dis) const override;

    /// get minimal node GID to be used in assign_degrees_of_freedom
    int get_minimal_node_gid_if_relevant(const Core::FE::Discretization& dis) const override;

    bool ignoreminnodegid_;  //< bool whether minnodegid is taken from the discretization or ignored

   private:
  };

}  // namespace Core::DOFSets

FOUR_C_NAMESPACE_CLOSE

#endif
