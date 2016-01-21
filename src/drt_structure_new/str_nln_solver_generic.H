/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_generic.H

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_GENERIC_H_
#define STR_NLN_SOLVER_GENERIC_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace NOX {
namespace Abstract {
  class Group;
} // namespace Abstract
} // namespace NOX
namespace STR {
namespace TIMINT {
  class Implicit;
  class BaseDataGlobalState;
  class BaseDataSDyn;
  class NoxInterface;
} // namespace TIMINT
namespace NLN {
namespace SOLVER {
  class Generic
  {
  public:
    //! constructor
    Generic();

    //! destructor
    virtual ~Generic() {};

    //! initialization
    virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn,
        const Teuchos::RCP<STR::TIMINT::NoxInterface> noxinterface);

    virtual void Setup() = 0;

    //! Reset internal storage before the nonlinear solution starts
    virtual void Reset() = 0;

    //! Solve the non-linear problem
    virtual INPAR::STR::ConvergenceStatus Solve() = 0;

    /*! returns the nox group for external and internal use
     *
     *  The nox group has to be initialized in one of the derived Setup() routines beforehand. */
    NOX::Abstract::Group& SolutionGroup();
    const NOX::Abstract::Group& GetSolutionGroup() const;

  protected:
    //! Returns true if Init() has been called
    inline const bool& IsInit() const
    { return isinit_; };

    //! Returns true if Setup() has been called
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Check if Init() and Setup() have been called
    void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }
    void CheckInit() const
    {
      if (!IsInit())
        dserror("You have to call Init() first!");
    }

    // Returns the global state data container pointer
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> DataGlobalStatePtr()
    {
      CheckInit();
      return gstate_ptr_;
    }

    // Returns the global state data container (read-only)
    const STR::TIMINT::BaseDataGlobalState& DataGlobalState() const
    {
      CheckInit();
      return *gstate_ptr_;
    }

    // Returns the global state data container (read and write)
    STR::TIMINT::BaseDataGlobalState& DataGlobalState()
    {
      CheckInit();
      return *gstate_ptr_;
    }

    // Returns the structural dynamics data container pointer
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> DataSDynPtr()
    {
      CheckInit();
      return sdyn_ptr_;
    }

    // Returns the structural dynamics data container (read-only)
    const STR::TIMINT::BaseDataSDyn& DataSDyn() const
    {
      CheckInit();
      return *sdyn_ptr_;
    }

    // Returns the structural dynamics data container (read and write)
    STR::TIMINT::BaseDataSDyn& DataSDyn()
    {
      CheckInit();
      return *sdyn_ptr_;
    }

    // Returns the non-linear solver implicit time integration interface pointer
    const Teuchos::RCP<STR::TIMINT::NoxInterface> NoxInterfacePtr()
    {
      CheckInit();
      return noxinterface_ptr_;
    }

    // Returns the non-linear solver implicit time integration interface (read-only)
    const STR::TIMINT::NoxInterface& NoxInterface() const
    {
      CheckInit();
      return *noxinterface_ptr_;
    }

    // Returns the non-linear solver implicit time integration interface (read and write)
    STR::TIMINT::NoxInterface& NoxInterface()
    {
      CheckInit();
      return *noxinterface_ptr_;
    }

    /*! returns the nox group (pointer) (only for internal use)
     *
     *  The nox group has to be initialized in one of the derived Setup() routines. */
    NOX::Abstract::Group& Group();
    Teuchos::RCP<NOX::Abstract::Group>& GroupPtr();

  protected:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;

  private:
    //! global state data container of the time integrator
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! structural dynamics data container of the time integrator
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn_ptr_;

    //! required interface pointer to the implicit time integrator (call back)
    Teuchos::RCP<STR::TIMINT::NoxInterface> noxinterface_ptr_;

    //! nox group
    Teuchos::RCP<NOX::Abstract::Group> group_ptr_;
  };  // class Generic
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_GENERIC_H_ */
