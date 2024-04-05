/*-----------------------------------------------------------*/
/*! \file
\brief factory for time integration base strategy and data container


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_FACTORY_HPP


#include "baci_config.hpp"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace STR
{
  class TimAda;
  namespace TIMINT
  {
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
      virtual ~Factory() = default;

      //! Build the implicit or explicit time integration strategies
      Teuchos::RCP<STR::TIMINT::Base> BuildStrategy(const Teuchos::ParameterList& sdyn) const;

      //! Build the structural dynamics data container
      Teuchos::RCP<STR::TIMINT::BaseDataSDyn> BuildDataSDyn(
          const Teuchos::ParameterList& sdyn) const;

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
        const Teuchos::ParameterList& ioflags,       //!< input-output-flags
        const Teuchos::ParameterList& sdyn,          //!< structural dynamic flags
        const Teuchos::ParameterList& xparams,       //!< extra flags
        const Teuchos::ParameterList& taflags,       //!< adaptive input flags
        Teuchos::RCP<STR::TIMINT::Base> ti_strategy  //!< marching time integrator
    );

    /*! Non-member function, which relates to the STR::TIMINT::Factory
     *  Please call this method, if you want to build a new structural dynamics data container. */
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> BuildDataSDyn(const Teuchos::ParameterList& sdyn);

    /*! Non-member function, which relates to the STR::TIMINT::Factory
     *  Please call this method, if you want to build a new global state data container. */
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> BuildDataGlobalState();
  }  // namespace TIMINT
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif
