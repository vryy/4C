/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate constraint enforcement implementations for Fluid-beam
interaction


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FBI_CONSTRAINTENFORCER_FACTORY_HPP
#define FOUR_C_FBI_CONSTRAINTENFORCER_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FBIConstraintenforcer;

  /**
   *  \brief Factory that creates the appropriate constraint enforcement strategy
   *
   *  To create a constraint enforcer implementation for FBI, call the static CreateEnforcer
   * function directly! No instance of ConstraintEnforcerFactory has to be created! If you try to
   * call the constructor, you will get an error message, since it is set to be private.
   */
  class ConstraintEnforcerFactory
  {
   private:
    /// constructor
    ConstraintEnforcerFactory() = delete;

   public:
    /**
     *  \brief Creates the appropriate constraint enforcement strategy
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] fsidyn List of FSI Input parameters
     * \param[in] fbidyn List of FBI Input parameters
     *
     * \return FBI constraint enforcement strategy
     */
    static Teuchos::RCP<FBIConstraintenforcer> create_enforcer(
        const Teuchos::ParameterList& fsidyn, const Teuchos::ParameterList& fbidyn);
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
