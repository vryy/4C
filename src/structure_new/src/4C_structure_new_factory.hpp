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

namespace Solid
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
    Teuchos::RCP<Solid::Integrator> build_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the desired  dirichlet boundary condition object
    Teuchos::RCP<Solid::Dbc> build_dbc(const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

   protected:
    //! build the implicit integrator
    Teuchos::RCP<Solid::Integrator> build_implicit_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the explicit integrator
    Teuchos::RCP<Solid::Integrator> build_explicit_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

  };  // class Factory

  /*! \brief Non-member function, which relates to the Solid::Factory
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<Solid::Integrator> build_integrator(const Solid::TimeInt::BaseDataSDyn& datasdyn);

  /*! \brief Non-member function, which relates to the Solid::Factory class
   *
   * \note Call this method from outside!
   */
  Teuchos::RCP<Solid::Dbc> build_dbc(const Solid::TimeInt::BaseDataSDyn& datasdyn);
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
