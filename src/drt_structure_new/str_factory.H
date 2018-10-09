/*-----------------------------------------------------------*/
/*!
\file str_factory.H

\brief factory for time integrator

\maintainer Michael Hiermeier

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_FACTORY_H_
#define SRC_DRT_STRUCTURE_NEW_STR_FACTORY_H_

#include <Teuchos_RCP.hpp>

// forward declarations
namespace ADAPTER
{
  class StructureBaseAlgorithmNew;
}  // namespace ADAPTER

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace STR
{
  class Integrator;
  class Dbc;

  namespace TIMINT
  {
    class BaseDataSDyn;
  }  // namespace TIMINT

  /*! \brief Factory to build the desired implicit/explicit integrator
   *
   *  \author Michael Hiermeier */
  class Factory
  {
   public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory(){};

    //! build the internal integrator
    Teuchos::RCP<STR::Integrator> BuildIntegrator(const STR::TIMINT::BaseDataSDyn& datasdyn) const;

    //! build the desired  dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> BuildDbc(const STR::TIMINT::BaseDataSDyn& datasdyn) const;

   protected:
    //! build the implicit integrator
    Teuchos::RCP<STR::Integrator> BuildImplicitIntegrator(
        const STR::TIMINT::BaseDataSDyn& datasdyn) const;

    //! build the explicit integrator
    Teuchos::RCP<STR::Integrator> BuildExplicitIntegrator(
        const STR::TIMINT::BaseDataSDyn& datasdyn) const;

  };  // class Factory

  /*! \brief Non-member function, which relates to the STR::Factory
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<STR::Integrator> BuildIntegrator(const STR::TIMINT::BaseDataSDyn& datasdyn);

  /*! \brief Non-member function, which relates to the STR::Factory class
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<STR::Dbc> BuildDbc(const STR::TIMINT::BaseDataSDyn& datasdyn);
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_FACTORY_H_ */
