/*-----------------------------------------------------------*/
/*!

\brief Generic class of the non-linear structural solvers.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_GENERIC_H_
#define STR_NLN_SOLVER_GENERIC_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX
namespace STR
{
  class Integrator;
  namespace TIMINT
  {
    class Implicit;
    class BaseDataGlobalState;
    class BaseDataSDyn;
    class Base;
    class NoxInterface;
  }  // namespace TIMINT
  namespace NLN
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
        virtual ~Generic(){};

        //! initialization
        virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate,
            const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn,
            const Teuchos::RCP<STR::TIMINT::NoxInterface>& noxinterface,
            const Teuchos::RCP<STR::Integrator>& integrator,
            const Teuchos::RCP<const STR::TIMINT::Base>& timint);

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
        virtual INPAR::STR::ConvergenceStatus Solve() = 0;

        /*! returns the nox group for external and internal use
         *
         *  The nox group has to be initialized in one of the derived Setup() routines beforehand.
         */
        NOX::Abstract::Group& SolutionGroup();
        const NOX::Abstract::Group& GetSolutionGroup() const;

        //! Get the number of nonlinear iterations
        virtual int GetNumNlnIterations() const = 0;

       protected:
        //! Returns true if Init() has been called
        inline const bool& IsInit() const { return isinit_; };

        //! Returns true if Setup() has been called
        inline const bool& IsSetup() const { return issetup_; };

        //! Check if Init() and Setup() have been called
        void CheckInitSetup() const
        {
          if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
        }

        //! Check if Init() has been called
        void CheckInit() const
        {
          if (!IsInit()) dserror("You have to call Init() first!");
        }

        //! Returns the global state data container pointer
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> DataGlobalStatePtr()
        {
          CheckInit();
          return gstate_ptr_;
        }

        //! Returns the global state data container (read-only)
        const STR::TIMINT::BaseDataGlobalState& DataGlobalState() const
        {
          CheckInit();
          return *gstate_ptr_;
        }

        //! Returns the global state data container (read and write)
        STR::TIMINT::BaseDataGlobalState& DataGlobalState()
        {
          CheckInit();
          return *gstate_ptr_;
        }

        //! Returns the structural dynamics data container pointer
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> DataSDynPtr()
        {
          CheckInit();
          return sdyn_ptr_;
        }

        //! Returns the structural dynamics data container (read-only)
        const STR::TIMINT::BaseDataSDyn& DataSDyn() const
        {
          CheckInit();
          return *sdyn_ptr_;
        }

        //! Returns the structural dynamics data container (read and write)
        STR::TIMINT::BaseDataSDyn& DataSDyn()
        {
          CheckInit();
          return *sdyn_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface pointer
        const Teuchos::RCP<STR::TIMINT::NoxInterface> NoxInterfacePtr()
        {
          CheckInit();
          return noxinterface_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface (read-only)
        const STR::TIMINT::NoxInterface& NoxInterface() const
        {
          CheckInit();
          return *noxinterface_ptr_;
        }

        //! Returns the non-linear solver implicit time integration interface (read and write)
        STR::TIMINT::NoxInterface& NoxInterface()
        {
          CheckInit();
          return *noxinterface_ptr_;
        }

        STR::Integrator& Integrator()
        {
          CheckInit();
          return *int_ptr_;
        }

        const STR::Integrator& Integrator() const
        {
          CheckInit();
          return *int_ptr_;
        }

        //! Returns the underlying time integration strategy
        const STR::TIMINT::Base& TimInt() const
        {
          CheckInit();
          return *timint_ptr_;
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

        //! pointer to the current time integrator
        Teuchos::RCP<STR::Integrator> int_ptr_;

        //! pointer to the time integration strategy
        Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

        //! nox group
        Teuchos::RCP<NOX::Abstract::Group> group_ptr_;

      };  // namespace SOLVER
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif /* STR_NLN_SOLVER_GENERIC_H_ */
