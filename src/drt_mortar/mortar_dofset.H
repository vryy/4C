/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of degrees of freedom special for mortar coupling

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef MORTAR_DOFSET_H
#define MORTAR_DOFSET_H

#include "../drt_lib/drt_dofset.H"


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
    virtual ~MortarDofSet() = default;

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new MortarDofSet(*this)); }

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
    virtual int AssignDegreesOfFreedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start);

    //! @}

   protected:
  };  // class MortarDofSet
}  // namespace MORTAR

#endif  // #ifndef MORTAR_DOFSET_H
