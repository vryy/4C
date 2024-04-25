/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate immersed coupling implementations for Fluid-beam
interaction


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FBI_IMMERSEDCOUPLER_FACTORY_HPP
#define FOUR_C_FBI_IMMERSEDCOUPLER_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace FBI
{
  class FBIGeometryCoupler;

  /**
   *  \brief Factory that creates the appropriate geometry coupler strategy
   *
   *  To create an immersed coupling implementation for FBI, call the static CreateGeometryCoupler
   * function directly! No instance of GeometryCouplerFactory has to be created! If you try to
   * call the constructor, you will get an error message.
   */
  class GeometryCouplerFactory
  {
    /// constructor
    GeometryCouplerFactory() = delete;

   public:
    /**
     *  \brief Creates the appropriate geometry coupler strategy
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] fbidyn List of FBI Input parameters
     *
     * \return FBI geometry coupler strategy
     */
    static Teuchos::RCP<FBI::FBIGeometryCoupler> CreateGeometryCoupler(
        const Teuchos::ParameterList& fbidyn);
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
