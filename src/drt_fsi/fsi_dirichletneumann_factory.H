/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate DirichletNeumann Algorithm

\maintainer Nora Hagmeyer

\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_FACTORY_H
#define FSI_DIRICHLETNEUMANN_FACTORY_H

#include "fsi_dirichletneumann.H"

namespace FSI
{
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
    /// destructor
    ~DirichletNeumannFactory(){};

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

#endif
