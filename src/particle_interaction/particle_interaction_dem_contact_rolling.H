/*---------------------------------------------------------------------------*/
/*! \file
\brief rolling contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H
#define PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactRollingBase
  {
   public:
    //! constructor
    explicit DEMContactRollingBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactRollingBase() = default;

    //! init rolling contact handler
    virtual void Init();

    //! setup rolling contact handler
    virtual void Setup(const double& k_normal);

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
        const double& mu_rolling, const double& normalcontactforce,
        double* rollingcontactmoment) const = 0;

    //! evaluate rolling potential energy
    virtual void RollingPotentialEnergy(
        const double* gap_rolling, double& rollingpotentialenergy) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! timestep size
    double dt_;

    //! coefficient of restitution
    const double e_;

    //! particle Poisson ratio
    const double nue_;

    //! rolling contact damping factor
    double d_rolling_fac_;
  };

  class DEMContactRollingViscous : public DEMContactRollingBase
  {
   public:
    //! constructor
    explicit DEMContactRollingViscous(const Teuchos::ParameterList& params);

    //! init rolling contact handler
    void Init() override;

    //! setup rolling contact handler
    void Setup(const double& k_normal) override;

    //! calculate effective radius
    void EffectiveRadiusParticle(const double* radius_i, const double* radius_j, const double& gap,
        double& r_eff) const override;

    //! calculate relative rolling velocity
    void RelativeRollingVelocity(const double& r_eff, const double* normal, const double* angvel_i,
        const double* angvel_j, double* v_rel_rolling) const override;

    //! calculate rolling contact moment
    void RollingContactMoment(double* gap_rolling, bool& stick_rolling, const double* normal,
        const double* v_rel_rolling, const double& m_eff, const double& r_eff,
        const double& mu_rolling, const double& normalcontactforce,
        double* rollingcontactmoment) const override;

    //! evaluate rolling potential energy
    void RollingPotentialEnergy(
        const double* gap_rolling, double& rollingpotentialenergy) const override;

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

    //! setup rolling contact handler
    void Setup(const double& k_normal) override;

    //! calculate effective radius
    void EffectiveRadiusParticle(const double* radius_i, const double* radius_j, const double& gap,
        double& r_eff) const override;

    //! calculate relative rolling velocity
    void RelativeRollingVelocity(const double& r_eff, const double* normal, const double* angvel_i,
        const double* angvel_j, double* v_rel_rolling) const override;

    //! calculate rolling contact moment
    void RollingContactMoment(double* gap_rolling, bool& stick_rolling, const double* normal,
        const double* v_rel_rolling, const double& m_eff, const double& r_eff,
        const double& mu_rolling, const double& normalcontactforce,
        double* rollingcontactmoment) const override;

    //! evaluate rolling potential energy
    void RollingPotentialEnergy(
        const double* gap_rolling, double& rollingpotentialenergy) const override;

   private:
    //! rolling contact stiffness
    double k_rolling_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
