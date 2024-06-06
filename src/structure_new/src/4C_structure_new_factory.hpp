/*-----------------------------------------------------------*/
/*! \file

\brief factory for time integrator


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class StructureBaseAlgorithmNew;
}  // namespace Adapter

namespace STR
{
  class Integrator;
  class Dbc;

  namespace TimeInt
  {
    class BaseDataSDyn;
  }  // namespace TimeInt

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
    Teuchos::RCP<STR::Integrator> build_integrator(
        const STR::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the desired  dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> build_dbc(const STR::TimeInt::BaseDataSDyn& datasdyn) const;

   protected:
    //! build the implicit integrator
    Teuchos::RCP<STR::Integrator> build_implicit_integrator(
        const STR::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the explicit integrator
    Teuchos::RCP<STR::Integrator> build_explicit_integrator(
        const STR::TimeInt::BaseDataSDyn& datasdyn) const;

  };  // class Factory

  /*! \brief Non-member function, which relates to the STR::Factory
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<STR::Integrator> build_integrator(const STR::TimeInt::BaseDataSDyn& datasdyn);

  /*! \brief Non-member function, which relates to the STR::Factory class
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<STR::Dbc> build_dbc(const STR::TimeInt::BaseDataSDyn& datasdyn);
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
