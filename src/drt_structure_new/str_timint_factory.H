/*
 * str_timint_factory.H
 *
 *  Created on: Sep 4, 2015
 *      Author: hiermeier
 */

#ifndef STR_TIMINT_FACTORY_H_
#define STR_TIMINT_FACTORY_H_


#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos {
  class ParameterList;
}
namespace STR {
namespace TIMINT {
  class Base;
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    //! Build the implicit or explicit time integrator
    Teuchos::RCP<STR::TIMINT::Base> BuildTimeIntegrator(
        const Teuchos::ParameterList& sdyn) const;

  private:
    //! Build the implicit time integrator
    Teuchos::RCP<STR::TIMINT::Base> BuildImplicitTimeIntegrator(
        const Teuchos::ParameterList& sdyn) const;

    //! Build the explicit time integrator
    Teuchos::RCP<STR::TIMINT::Base> BuildExplicitTimeIntegrator(
        const Teuchos::ParameterList& sdyn) const;
  };

  //! non-member function
  Teuchos::RCP<STR::TIMINT::Base> BuildTimeIntegrator(const Teuchos::ParameterList& sdyn);
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_FACTORY_H_ */
