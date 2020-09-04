/*---------------------------------------------------------------------------*/
/*! \file
\brief dirichlet boundary condition handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_DIRICHLET_BC_H
#define PARTICLE_DIRICHLET_BC_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief dirichlet boundary condition handler for particle simulations
   *
   * \author Sebastian Fuchs \date 07/2018
   */
  class DirichletBoundaryConditionHandler
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] params particle simulation parameter list
     */
    explicit DirichletBoundaryConditionHandler(const Teuchos::ParameterList& params);

    //! destructor
    ~DirichletBoundaryConditionHandler() = default;

    /*!
     * \brief init dirichlet boundary condition handler
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void Init();

    /*!
     * \brief setup dirichlet boundary condition handler
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief get reference to set of particle types subjected to dirichlet boundary conditions
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \return set of particle types subjected to dirichlet boundary conditions
     */
    const std::set<PARTICLEENGINE::TypeEnum>& GetParticleTypesSubjectedToDirichletBCSet() const
    {
      return typessubjectedtodirichletbc_;
    };

    /*!
     * \brief insert dirichlet boundary condition dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    /*!
     * \brief set particle reference position
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void SetParticleReferencePosition() const;

    /*!
     * \brief evaluate dirichlet boundary condition
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] evaltime evaluation time
     * \param[in] evalpos  flag to indicate evaluation of position
     * \param[in] evalvel  flag to indicate evaluation of velocity
     * \param[in] evalacc  flag to indicate evaluation of acceleration
     */
    void EvaluateDirichletBoundaryCondition(
        const double& evaltime, const bool evalpos, const bool evalvel, const bool evalacc) const;

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! relating particle types to function ids of dirichlet boundary conditions
    std::map<PARTICLEENGINE::TypeEnum, int> dirichletbctypetofunctid_;

    //! set of particle types subjected to dirichlet boundary conditions
    std::set<PARTICLEENGINE::TypeEnum> typessubjectedtodirichletbc_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
