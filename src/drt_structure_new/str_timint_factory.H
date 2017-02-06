/*-----------------------------------------------------------*/
/*!
\file str_timint_factory.H

\maintainer Michael Hiermeier

\date Sep 4, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_TIMINT_FACTORY_H_
#define STR_TIMINT_FACTORY_H_


#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos {
  class ParameterList;
}
namespace STR {
  class TimAda;
namespace TIMINT {
  class Base;
  class BaseDataSDyn;
  class BaseDataGlobalState;

  /*! \brief Factory to build the desired time integration strategy and
   *  the desired adaptive wrapper object
   *
   *  \author Michael Hiermeier */
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    //! Build the implicit or explicit time integration strategies
    Teuchos::RCP<STR::TIMINT::Base> BuildStrategy(
        const Teuchos::ParameterList& sdyn) const;

    //! Build adaptive auxiliary time integration strategy
    Teuchos::RCP<STR::TimAda> BuildAdaptiveWrapper(
        const Teuchos::ParameterList& ioflags,         //!< input-output-flags
        const Teuchos::ParameterList& sdyn,            //!< structural dynamic flags
        const Teuchos::ParameterList& xparams,         //!< extra flags
        const Teuchos::ParameterList& taflags,         //!< adaptive input flags
        Teuchos::RCP<STR::TIMINT::Base> ti_strategy    //!< marching time integrator
        ) const;

    //! Build the structural dynamics data container
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> BuildDataSDyn(
        const Teuchos::ParameterList& sdyn
        ) const;

    //! Build the global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> BuildDataGlobalState() const;

  private:
    //! Build the implicit time integration strategies
    Teuchos::RCP<STR::TIMINT::Base> BuildImplicitStrategy(
        const Teuchos::ParameterList& sdyn) const;

    //! Build the explicit time integration strategies
    Teuchos::RCP<STR::TIMINT::Base> BuildExplicitStrategy(
        const Teuchos::ParameterList& sdyn) const;
  };  // class Factory

  /*! Non-member function, which relates to the STR::TIMINT::Factory
   *  Please call this method, if you want to build a new time integration strategy. */
  Teuchos::RCP<STR::TIMINT::Base> BuildStrategy(const Teuchos::ParameterList& sdyn);

  /*! Non-member function, which relates to the STR::TIMINT::Factory
   *  Please call this method, if you want to build a new adaptive wrapper object. */
  Teuchos::RCP<STR::TimAda> BuildAdaptiveWrapper(
      const Teuchos::ParameterList& ioflags,         //!< input-output-flags
      const Teuchos::ParameterList& sdyn,            //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,         //!< extra flags
      const Teuchos::ParameterList& taflags,         //!< adaptive input flags
      Teuchos::RCP<STR::TIMINT::Base> ti_strategy    //!< marching time integrator
      );

  /*! Non-member function, which relates to the STR::TIMINT::Factory
   *  Please call this method, if you want to build a new structural dynamics data container. */
  Teuchos::RCP<STR::TIMINT::BaseDataSDyn> BuildDataSDyn(
      const Teuchos::ParameterList& sdyn
      );

  /*! Non-member function, which relates to the STR::TIMINT::Factory
   *  Please call this method, if you want to build a new global state data container. */
  Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> BuildDataGlobalState();
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_FACTORY_H_ */
