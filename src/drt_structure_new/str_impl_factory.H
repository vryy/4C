/*
 * str_impl_factory.H
 *
 *  Created on: Sep 2, 2015
 *      Author: hiermeier
 */

#ifndef STR_IMPL_FACTORY_H_
#define STR_IMPL_FACTORY_H_

#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace STR {
namespace IMPLICIT {
  class Generic;

  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    //! Build specified factory
    Teuchos::RCP<STR::IMPLICIT::Generic> BuildImplicitIntegrator(
        const enum INPAR::STR::DynamicType& dynType,
        const enum INPAR::STR::PreStress& preStressType) const;
  };

  // non-member function
  Teuchos::RCP<STR::IMPLICIT::Generic> BuildImplicitIntegrator(
      const enum INPAR::STR::DynamicType& intType,
      const enum INPAR::STR::PreStress& preStressType);
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_FACTORY_H_ */
