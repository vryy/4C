/*-----------------------------------------------------------*/
/*! \file

\brief factory for time integrator


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_FACTORY_H
#define BACI_STRUCTURE_NEW_FACTORY_H

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
    virtual ~Factory() = default;

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


#endif  // STRUCTURE_NEW_FACTORY_H
