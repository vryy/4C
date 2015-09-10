/*
 * str_nln_solver_generic.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_GENERIC_H_
#define STR_NLN_SOLVER_GENERIC_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace STR {
namespace NLN {
namespace SOLVER {
namespace INTERFACE {
  class Required;
} // namespace INTERFACE

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
        const Teuchos::RCP<STR::NLN::SOLVER::INTERFACE::Required> ireq);

    virtual void Setup() = 0;

    //! Reset internal storage before the nonlinear solution starts
    virtual void Reset() = 0;

    //! Solve the non-linear problem
    virtual void Solve() = 0;

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

    //! required interface pointer
    Teuchos::RCP<STR::NLN::SOLVER::INTERFACE::Required> ireq_;
  };
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_GENERIC_H_ */
