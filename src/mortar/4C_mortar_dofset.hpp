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
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  /*!
  \brief A set of degrees of freedom special for mortar coupling

  */
  class DofSet : public CORE::Dofsets::DofSet
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    /*!
    \brief Standard Constructor

    */
    DofSet();



    /// create a copy of this object
    Teuchos::RCP<CORE::Dofsets::DofSet> Clone() override
    {
      return Teuchos::rcp(new MORTAR::DofSet(*this));
    }

    //! @}

    //! @name Construction
    //! @{

    /*!
    \brief Assign dof numbers to all elements and nodes of the discretization

    This specialized DofSet does not generate DOFs by itself, but rather extracts them from the
    mortar nodes and sticks them into the discretization of the mortar interface.

    @param[in] dis Discretization of a mortar interface
    @param[in] dspos Position of DOfSet inside its discretization
    @param[in] start User-defined offset for DOF numbering [currently not supported]

    @return Maximum dof number of this dofset
    */
    int assign_degrees_of_freedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start) override;

    //! @}

   protected:
  };  // class DofSet
}  // namespace MORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
