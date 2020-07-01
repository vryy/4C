/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate constraint enforcement implementations for Fluid-beam
interaction


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FBI_CONSTRAINTENFORCER_FACTORY_H
#define FSI_CONSTRAINTENFORCER_FACTORY_H

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}
namespace ADAPTER
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
    ConstraintEnforcerFactory(){};

   public:
    /// destructor
    ~ConstraintEnforcerFactory(){};

    /**
     *  \brief Creates the appropriate constraint enforcement strategy
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] fsidyn List of FSI Input parameters
     *
     * \return FBI constraint enforcement strategy
     */
    static Teuchos::RCP<FBIConstraintenforcer> CreateEnforcer(const Teuchos::ParameterList &fsidyn);
  };
}  // namespace ADAPTER

#endif
