/*---------------------------------------------------------------------*/
/*! \file

\brief This class is inherited from dofset and replaces the
       method assign_degrees_of_freedom by a version that does not query
       the static_dofsets_ list for the max GID, but always starts from
       0. Also, it does not register with the static_dofsets_ list. This
       class is intended to be used for xfem approaches. It provides a
       dofset without xfem dofs for output routines.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_DOFSET_INDEPENDENT_HPP
#define FOUR_C_FEM_DOFSET_INDEPENDENT_HPP

#include "4C_config.hpp"

#include "4C_fem_dofset.hpp"

#include <Teuchos_RCP.hpp>

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
    Teuchos::RCP<DofSet> Clone() override { return Teuchos::rcp(new IndependentDofSet(*this)); }

    /// Add Dof Set to list #static_dofsets_
    void AddDofSettoList() override;

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
