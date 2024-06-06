/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of degrees of freedom special for mortar coupling

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_DOFSET_HPP
#define FOUR_C_MORTAR_DOFSET_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;
}

namespace Mortar
{
  /*!
  \brief A set of degrees of freedom special for mortar coupling

  */
  class DofSet : public Core::DOFSets::DofSet
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    /*!
    \brief Standard Constructor

    */
    DofSet();



    /// create a copy of this object
    Teuchos::RCP<Core::DOFSets::DofSet> Clone() override
    {
      return Teuchos::rcp(new Mortar::DofSet(*this));
    }

    //! @}

    //! @name Construction
    //! @{

    /*!
    \brief Assign dof numbers to all elements and nodes of the discretization

    This specialized DofSet does not generate DOFs by itself, but rather extracts them from the
    mortar nodes and sticks them into the discretization of the mortar interface.

    @param[in] dis discretization of a mortar interface
    @param[in] dspos Position of DOfSet inside its discretization
    @param[in] start User-defined offset for DOF numbering [currently not supported]

    @return Maximum dof number of this dofset
    */
    int assign_degrees_of_freedom(
        const Discret::Discretization& dis, const unsigned dspos, const int start) override;

    //! @}

   protected:
  };  // class DofSet
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif
