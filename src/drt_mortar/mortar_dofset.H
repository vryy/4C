/*!----------------------------------------------------------------------
\file mortar_dofset.H

\brief A set of degrees of freedom special for mortar coupling

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_DOFSET_H
#define MORTAR_DOFSET_H

#include "../drt_lib/drt_dofset.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{
  /*!
  \brief A set of degrees of freedom special for mortar coupling

  */
  class MortarDofSet : public DRT::DofSet
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    */
    MortarDofSet();

    /*!
    \brief Destructor

    */
    virtual ~MortarDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new MortarDofSet(*this)); }

    //@}

    //! @name Construction

    /// Assign dof numbers to all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start);

    //@}

   protected:
  };  // class MortarDofSet
}  // namespace MORTAR

#endif  // #ifndef MORTAR_DOFSET_H
