/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of degrees of freedom special for mortar coupling

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef BACI_MORTAR_DOFSET_H
#define BACI_MORTAR_DOFSET_H

#include "baci_lib_dofset.H"


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
  class MortarDofSet : public DRT::DofSet
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    /*!
    \brief Standard Constructor

    */
    MortarDofSet();

    /*!
    \brief Destructor

    */
    ~MortarDofSet() override = default;

    /// create a copy of this object
    Teuchos::RCP<DofSet> Clone() override { return Teuchos::rcp(new MortarDofSet(*this)); }

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
    int AssignDegreesOfFreedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start) override;

    //! @}

   protected:
  };  // class MortarDofSet
}  // namespace MORTAR

#endif  // MORTAR_DOFSET_H
