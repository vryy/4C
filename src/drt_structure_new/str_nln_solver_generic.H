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
namespace STR {
namespace TIMINT {
  class Implicit;
  class BaseDataGlobalState;
  class BaseDataSDyn;
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
        const Teuchos::RCP<STR::TIMINT::Implicit> implicit);

    virtual void Setup() = 0;

    //! Reset internal storage before the nonlinear solution starts
    virtual void Reset() = 0;

    //! Solve the non-linear problem
    virtual INPAR::STR::ConvergenceStatus Solve() = 0;

  protected:
    //! Returns true if Init() has been called
    const bool& IsInit() const
    { return isinit_; };

    //! Returns true if Setup() has been called
    const bool& IsSetup() const
    { return issetup_; };

    //! Check if Init() and Setup() have been called
    virtual void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }
    virtual void CheckInit() const
    {
      if (!IsInit())
        dserror("You have to call Init() first!");
    }

    // Returns the global state data container pointer
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> DataGlobalStatePtr()
    {
      CheckInitSetup();
      return gstate_;
    }

    // Returns the global state data container (read-only)
    const STR::TIMINT::BaseDataGlobalState& DataGlobalState() const
    {
      CheckInitSetup();
      return *gstate_;
    }

    // Returns the global state data container (read and write)
    STR::TIMINT::BaseDataGlobalState& DataGlobalState()
    {
      CheckInitSetup();
      return *gstate_;
    }

    // Returns the structural dynamics data container pointer
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> DataSDynPtr()
    {
      CheckInitSetup();
      return sdyn_;
    }

    // Returns the structural dynamics data container (read-only)
    const STR::TIMINT::BaseDataSDyn& DataSDyn() const
    {
      CheckInitSetup();
      return *sdyn_;
    }

    // Returns the structural dynamics data container (read and write)
    STR::TIMINT::BaseDataSDyn& DataSDyn()
    {
      CheckInitSetup();
      return *sdyn_;
    }

    // Returns the non-linear solver implicit time integration interface pointer
    const Teuchos::RCP<STR::TIMINT::Implicit> ImplicitTimIntPtr()
    {
      CheckInitSetup();
      return iimpl_;
    }

    // Returns the non-linear solver implicit time integration interface (read-only)
    const STR::TIMINT::Implicit& ImplicitTimInt() const
    {
      CheckInitSetup();
      return *iimpl_;
    }

    // Returns the non-linear solver implicit time integration interface (read and write)
    STR::TIMINT::Implicit& ImplicitTimInt()
    {
      CheckInitSetup();
      return *iimpl_;
    }

  protected:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;

  private:
    //! global state data container of the time integrator
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_;

    //! structural dynamics data container of the time integrator
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn_;

    //! required interface pointer to the implicit time integrator (call back)
    Teuchos::RCP<STR::TIMINT::Implicit> iimpl_;
  };  // class Generic
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_GENERIC_H_ */
