/*---------------------------------------------------------------------------*/
/*!
\brief rolling contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H
#define PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactRollingBase
  {
   public:
    //! constructor
    explicit DEMContactRollingBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactRollingBase(){};

    //! init rolling contact handler
    virtual void Init();

    //! setup rolling contact handler
    virtual void Setup(const double& k_normal);

    //! write restart of rolling contact handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of rolling contact handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize) final;

    //! calculate effective radius
    virtual void EffectiveRadiusParticle(
        const double* radius_i, const double* radius_j, const double& gap, double& r_eff) const = 0;

    //! calculate relative rolling velocity
    virtual void RelativeRollingVelocity(const double& r_eff, const double* normal,
        const double* angvel_i, const double* angvel_j, double* vel_rel_rolling) const = 0;

    //! calculate rolling contact moment
    virtual void RollingContactMoment(double* gap_rolling, bool& stick_rolling,
        const double* normal, const double* v_rel_rolling, const double& m_eff, const double& r_eff,
        const double& normalcontactforce, double* rollingcontactmoment) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! timestep size
    double dt_;

    //! coefficient of restitution (particle-particle)
    const double e_;

    //! particle Poisson ratio
    const double nue_;

    //! rolling contact friction coefficient (particle-particle)
    const double mu_rolling_;

    //! rolling contact damping factor (particle-particle)
    double d_rolling_fac_;
  };

  class DEMContactRollingViscous : public DEMContactRollingBase
  {
   public:
    //! constructor
    explicit DEMContactRollingViscous(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactRollingViscous(){};

    //! setup rolling contact handler
    virtual void Setup(const double& k_normal);

    //! calculate effective radius
    virtual void EffectiveRadiusParticle(const double* radius_i, const double* radius_j,
        const double& gap, double& r_eff) const override;

    //! calculate relative rolling velocity
    virtual void RelativeRollingVelocity(const double& r_eff, const double* normal,
        const double* angvel_i, const double* angvel_j, double* v_rel_rolling) const override;

    //! calculate rolling contact moment
    virtual void RollingContactMoment(double* gap_rolling, bool& stick_rolling,
        const double* normal, const double* v_rel_rolling, const double& m_eff, const double& r_eff,
        const double& normalcontactforce, double* rollingcontactmoment) const override;

   private:
    //! particle Young's modulus
    const double young_;

    //! maximum expected particle velocity
    const double v_max_;
  };

  class DEMContactRollingCoulomb : public DEMContactRollingBase
  {
   public:
    //! constructor
    explicit DEMContactRollingCoulomb(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactRollingCoulomb(){};

    //! setup rolling contact handler
    virtual void Setup(const double& k_normal);

    //! calculate effective radius
    virtual void EffectiveRadiusParticle(const double* radius_i, const double* radius_j,
        const double& gap, double& r_eff) const override;

    //! calculate relative rolling velocity
    virtual void RelativeRollingVelocity(const double& r_eff, const double* normal,
        const double* angvel_i, const double* angvel_j, double* v_rel_rolling) const override;

    //! calculate rolling contact moment
    virtual void RollingContactMoment(double* gap_rolling, bool& stick_rolling,
        const double* normal, const double* v_rel_rolling, const double& m_eff, const double& r_eff,
        const double& normalcontactforce, double* rollingcontactmoment) const override;

   private:
    //! rolling contact stiffness (particle-particle)
    double k_rolling_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
