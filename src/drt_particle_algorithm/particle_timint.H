/*---------------------------------------------------------------------------*/
/*! \file
\brief time integration for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_H
#define PARTICLE_TIMINT_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEALGORITHM
{
  class DirichletBoundaryConditionHandler;
  class TemperatureBoundaryConditionHandler;
}  // namespace PARTICLEALGORITHM

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class TimInt
  {
   public:
    //! constructor
    explicit TimInt(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~TimInt();

    //! init particle time integration
    virtual void Init();

    //! setup particle time integration
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const bool modifiedstates);

    //! write restart of particle time integration
    void WriteRestart(const int step, const double time) const;

    //! read restart of particle time integration
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert integration dependent states of all particle types
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! time integration scheme specific initialization routine
    virtual void SetInitialStates();

    //! set current time
    virtual void SetCurrentTime(const double currenttime) final;

    //! time integration scheme specific pre-interaction routine
    virtual void PreInteractionRoutine() = 0;

    //! time integration scheme specific post-interaction routine
    virtual void PostInteractionRoutine() = 0;

   private:
    //! init dirichlet boundary condition handler
    void InitDirichletBoundaryCondition();

    //! init temperature boundary condition handler
    void InitTemperatureBoundaryCondition();

    //! add initial random noise to particle position
    void AddInitialRandomNoiseToPosition();

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! dirichlet boundary condition handler
    std::unique_ptr<PARTICLEALGORITHM::DirichletBoundaryConditionHandler>
        dirichletboundarycondition_;

    //! temperature boundary condition handler
    std::unique_ptr<PARTICLEALGORITHM::TemperatureBoundaryConditionHandler>
        temperatureboundarycondition_;

    //! set of particle types to integrate in time
    std::set<PARTICLEENGINE::TypeEnum> typestointegrate_;

    //! current time
    double time_;

    //! time step size
    double dt_;

    //! modified velocity and acceleration states
    bool modifiedstates_;
  };

  class TimIntSemiImplicitEuler : public TimInt
  {
    /*
     * semi-implicit Euler time integration scheme of first order accuracy (also denoted as
     * semi-explicit Euler or symplectic Euler scheme)
     *
     * v_{n+1} = v_{n} + dt * a_{n}      with    a_{n} = a( r_{n-1}, v_{n-1} )
     * x_{n+1} = x_{n} + dt * v_{n+1}
     *
     */

   public:
    //! constructor
    TimIntSemiImplicitEuler(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~TimIntSemiImplicitEuler(){};

    //! setup particle time integration
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const bool modifiedstates);

    //! time integration scheme specific pre-interaction routine
    virtual void PreInteractionRoutine() override;

    //! time integration scheme specific post-interaction routine
    virtual void PostInteractionRoutine() override;
  };

  class TimIntVelocityVerlet : public TimInt
  {
    /*
     * explicit velocity Verlet time integration scheme of second order accuracy (also denoted as
     * leapfrog scheme in kick-drift-kick form)
     *
     * v_{n+1/2} = v_{n}     + dt/2 * a_{n}        with    a_{n}   = a( r_{n}, v_{n-1/2} )
     * x_{n+1}   = x_{n}     + dt   * v_{n+1/2}
     * v_{n+1}   = v_{n+1/2} + dt/2 * a_{n+1}      with    a_{n+1} = a( r_{n+1}, v_{n+1/2} )
     *
     */

   public:
    //! constructor
    TimIntVelocityVerlet(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~TimIntVelocityVerlet(){};

    //! time integration scheme specific initialization routine
    virtual void SetInitialStates();

    //! time integration scheme specific pre-interaction routine
    virtual void PreInteractionRoutine() override;

    //! time integration scheme specific post-interaction routine
    virtual void PostInteractionRoutine() override;

   private:
    //! half time step size
    double dthalf_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
