/*-----------------------------------------------------------*/
/*! \file

\brief Generic class of the non-linear structural solvers.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_GENERIC_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_GENERIC_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  class Integrator;
  namespace TimeInt
  {
    class Implicit;
    class BaseDataGlobalState;
    class BaseDataSDyn;
    class Base;
    class NoxInterface;
  }  // namespace TimeInt
  namespace Nln
  {
    namespace SOLVER
    {
      /*! \brief Base class of all nonlinear solvers for structural dynamcis
       *
       */
      class Generic
      {
       public:
        //! constructor
        Generic();

        //! destructor
        virtual ~Generic() = default;

        //! initialization
        virtual void Init(const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate,
            const Teuchos::RCP<STR::TimeInt::BaseDataSDyn>& sdyn,
            const Teuchos::RCP<STR::TimeInt::NoxInterface>& noxinterface,
            const Teuchos::RCP<STR::Integrator>& integrator,
            const Teuchos::RCP<const STR::TimeInt::Base>& timint);

        //! Setup the nonlinear solver configuration
        virtual void Setup() = 0;

        /*! \brief Reset internal storage before the nonlinear solution starts
         *
         *  We actually (re-)build the nonlinear solver object here.
         *
         *  \warning It is not fully clear how rebuilding the nonlinear solver affects a possible
         *           re-use of the preconditioner for the linear system.
         */
        virtual void Reset() = 0;

        //! Solve the non-linear problem
        virtual Inpar::STR::ConvergenceStatus Solve() = 0;

        /*! returns the nox group for external and internal use
         *
         *  The nox group has to be initialized in one of the derived Setup() routines beforehand.
         */
        ::NOX::Abstract::Group& SolutionGroup();
        const ::NOX::Abstract::Group& get_solution_group() const;

        //! Get the number of nonlinear iterations
        virtual int GetNumNlnIterations() const = 0;

       protected:
        //! Returns true if Init() has been called
        inline const bool& is_init() const { return isinit_; };

        //! Returns true if Setup() has been called
        inline const bool& is_setup() const { return issetup_; };

        //! Check if Init() and Setup() have been called
        void check_init_setup() const
        {
          FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
        }

        //! Check if Init() has been called
        void check_init() const { FOUR_C_ASSERT(is_init(), "You have to call Init() first!"); }

        //! Returns the global state data container pointer
        Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> data_global_state_ptr()
        {
          check_init();
          return gstate_ptr_;
        }

        //! Returns the global state data container (read-only)
        const STR::TimeInt::BaseDataGlobalState& data_global_state() const
        {
          check_init();
          return *gstate_ptr_;
        }

        //! Returns the global state data container (read and write)
        STR::TimeInt::BaseDataGlobalState& data_global_state()
        {
          check_init();
          return *gstate_ptr_;
        }

        //! Returns the structural dynamics data container pointer
        Teuchos::RCP<STR::TimeInt::BaseDataSDyn> data_s_dyn_ptr()
        {
          check_init();
          return sdyn_ptr_;
        }

        //! Returns the structural dynamics data container (read-only)
        const STR::TimeInt::BaseDataSDyn& data_s_dyn() const
        {
          check_init();
          return *sdyn_ptr_;
        }

        //! Returns the structural dynamics data container (read and write)
        STR::TimeInt::BaseDataSDyn& data_s_dyn()
        {
          check_init();
          return *sdyn_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface pointer
        Teuchos::RCP<STR::TimeInt::NoxInterface> nox_interface_ptr()
        {
          check_init();
          return noxinterface_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface (read-only)
        const STR::TimeInt::NoxInterface& nox_interface() const
        {
          check_init();
          return *noxinterface_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface (read and write)
        STR::TimeInt::NoxInterface& nox_interface()
        {
          check_init();
          return *noxinterface_ptr_;
        }

        STR::Integrator& integrator()
        {
          check_init();
          return *int_ptr_;
        }

        const STR::Integrator& integrator() const
        {
          check_init();
          return *int_ptr_;
        }

        //! Returns the underlying time integration strategy
        const STR::TimeInt::Base& tim_int() const
        {
          check_init();
          return *timint_ptr_;
        }

        /*! returns the nox group (pointer) (only for internal use)
         *
         *  The nox group has to be initialized in one of the derived Setup() routines. */
        ::NOX::Abstract::Group& group();
        Teuchos::RCP<::NOX::Abstract::Group>& group_ptr();

       protected:
        //! init flag
        bool isinit_;

        //! setup flag
        bool issetup_;

       private:
        //! global state data container of the time integrator
        Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

        //! structural dynamics data container of the time integrator
        Teuchos::RCP<STR::TimeInt::BaseDataSDyn> sdyn_ptr_;

        //! required interface pointer to the implicit time integrator (call back)
        Teuchos::RCP<STR::TimeInt::NoxInterface> noxinterface_ptr_;

        //! pointer to the current time integrator
        Teuchos::RCP<STR::Integrator> int_ptr_;

        //! pointer to the time integration strategy
        Teuchos::RCP<const STR::TimeInt::Base> timint_ptr_;

        //! nox group
        Teuchos::RCP<::NOX::Abstract::Group> group_ptr_;

      };  // namespace SOLVER
    }     // namespace SOLVER
  }       // namespace Nln
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
