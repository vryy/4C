/*---------------------------------------------------------------------*/
/*! \file

\brief This class is inherited from dofset and replaces the
       method AssignDegreesOfFreedom by a version that does not query
       the static_dofsets_ list for the max GID, but always starts from
       0. Also, it does not register with the static_dofsets_ list. This
       class is intended to be used for xfem approaches. It provides a
       dofset without xfem dofs for output routines.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_DOFSET_INDEPENDENT_HPP
#define FOUR_C_DISCRETIZATION_DOFSET_INDEPENDENT_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CORE::Dofsets
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
    /// get first number to be used as Dof GID in AssignDegreesOfFreedom
    int GetFirstGIDNumberToBeUsed(const DRT::Discretization& dis) const override;

    /// get minimal node GID to be used in AssignDegreesOfFreedom
    int GetMinimalNodeGIDIfRelevant(const DRT::Discretization& dis) const override;

    bool ignoreminnodegid_;  //< bool whether minnodegid is taken from the discretization or ignored

   private:
  };

}  // namespace CORE::Dofsets

FOUR_C_NAMESPACE_CLOSE

#endif
