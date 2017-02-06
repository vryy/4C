/*---------------------------------------------------------------------------*/
/*!
\file xfield_state.H

\brief state management for one XFEM discretization

\level 3

\maintainer Michael Hiermeier

\date Jun 27, 2016
*/
/*---------------------------------------------------------------------------*/

#ifndef SRC_DRT_XFEM_XFIELD_STATE_H_
#define SRC_DRT_XFEM_XFIELD_STATE_H_

#include "../drt_lib/drt_dserror.H"
#include <Teuchos_RCP.hpp>

class Epetra_Map;

namespace DRT {
  class DiscretizationInterface;
} // namespace DRT
namespace GEO {
  class CutWizard;
} // namespace GEO

namespace XFEM {
  class ConditionManager;
  class XFEMDofSet;

  class XFieldState
  {
  public:
    /// constructor
    XFieldState();

    /// destructor
    virtual ~XFieldState() {};

    /** \brief initialize member variables for xfield<-->field couplings
     *
     *  Examples are the XFluidFluid problem and the XContact problem. */
    void Init(
        const Teuchos::RCP<XFEM::ConditionManager>&  condition_manager,
        const Teuchos::RCP<GEO::CutWizard>&          wizard,
        const Teuchos::RCP<XFEM::XFEMDofSet>&        xdofset,
        const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
        const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret);

    /// setup the stored state objects
    virtual void Setup() = 0;

    /// destroy the stored objects
    virtual bool Destroy() = 0;

    /// @name Accessors
    /// @{
    /// Get cut wizard
    GEO::CutWizard& CutWizard()
    {
      CheckInit();
      if (wizard_.is_null())
        dserror( "The CutWizard was not initialized! (Teuchos::null)" );
      return *wizard_;
    }

    /// Get condition manager
    XFEM::ConditionManager& ConditionManager()
    {
      CheckInit();
      if (condition_manager_.is_null())
        dserror( "The condition_manager was not initialized! (Teuchos::null)" );
      return *condition_manager_;
    }

    /// Get dofset of the cut discretization
    XFEM::XFEMDofSet& XDofSet()
    {
      CheckInit();
      if (xdofset_.is_null())
        dserror( "The xDoF set was not initialized! (Teuchos::null)" );
      return *xdofset_;
    }

  protected:
    /// Get cut wizard pointer
    Teuchos::RCP<GEO::CutWizard>& CutWizardPtr()
    {
      return wizard_;
    }

    /// Get condition manager pointer
    Teuchos::RCP<XFEM::ConditionManager>& ConditionManagerPtr()
    {
      return condition_manager_;
    }

    /// Get pointer to the dofset of the cut discretization
    Teuchos::RCP<XFEM::XFEMDofSet>& XDofSetPtr()
    {
      return xdofset_;
    }

    /// Returns the xFEM field discretizaton
    DRT::DiscretizationInterface& XFieldDiscret()
    {
      if (xfield_discret_ptr_.is_null())
        dserror("xfield_discret_ptr_ is NULL!");

      return *xfield_discret_ptr_;
    }

    /// Returns a pointer to the xFEM discretization
    Teuchos::RCP<DRT::DiscretizationInterface>& XFieldDiscretPtr()
    {
      return xfield_discret_ptr_;
    }

    /// Returns the standard field discretizaton
    DRT::DiscretizationInterface& FieldDiscret()
    {
      if (field_discret_ptr_.is_null())
        dserror("field_discret_ptr_ is NULL!");

      return *field_discret_ptr_;
    }

    /// Returns a pointer to the standard discretization
    Teuchos::RCP<DRT::DiscretizationInterface>& FieldDiscretPtr()
    {
      return field_discret_ptr_;
    }

    /// @}


  protected:
    //! check the initialization indicator
    inline void CheckInit() const
    {
      if (not isinit_)
        dserror("Call XFEM::XFieldState::Init() first!");
    }

    //! check the initialization and setup indicators
    inline void CheckInitSetup() const
    {
      if (not issetup_ or not isinit_)
        dserror("Call Init() and Setup() first!");
    }

  protected:
    /// init indicator
    bool isinit_;

    /// setup indicator
    bool issetup_;

  private:
    /// cut wizard
    Teuchos::RCP<GEO::CutWizard> wizard_;

    /// condition manager
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    /// XFEM dofset
    Teuchos::RCP<XFEM::XFEMDofSet> xdofset_;

    /// XFEM field discretization pointer
    Teuchos::RCP<DRT::DiscretizationInterface> xfield_discret_ptr_;

    /// field discretization pointer
    Teuchos::RCP<DRT::DiscretizationInterface> field_discret_ptr_;
  };  // class XFieldState
} // namespace XFEM


#endif /* SRC_DRT_XFEM_XFIELD_STATE_H_ */
