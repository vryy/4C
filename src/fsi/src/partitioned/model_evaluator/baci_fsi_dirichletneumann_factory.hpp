/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate DirichletNeumann Algorithm


\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_DIRICHLETNEUMANN_FACTORY_HPP
#define FOUR_C_FSI_DIRICHLETNEUMANN_FACTORY_HPP

#include "baci_config.hpp"

#include <Teuchos_RCPDecl.hpp>

class Epetra_Comm;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

namespace FSI
{
  class DirichletNeumann;
  /**
   *  \brief Factory that creates the appropriate DirichletNeumann algorithm
   *
   *  To create a DirichletNeumann algorithm, call the static CreateAlgorithm function directly! No
   * instance of DirichletNeumannFactory has to be created! If you try to call the constructor, you
   * will get an error message, since it is set to be private.
   */
  class DirichletNeumannFactory
  {
   private:
    /// constructor
    DirichletNeumannFactory(){};

   public:
    /**
     *  \brief Creates the appropriate DirichletNeumann algorithm
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] comm Epetra Communicator used in FSI::Partitioned for Terminal Output
     * \param[in] fsidyn List of FSI Input parameters
     *
     * \return Coupling algorithm based on Dirichlet-Neumann partitioning
     */
    static Teuchos::RCP<DirichletNeumann> CreateAlgorithm(
        const Epetra_Comm &comm, const Teuchos::ParameterList &fsidyn);
  };
}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif
