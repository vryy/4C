/*---------------------------------------------------------------------------*/
/*! \file
\brief tangential contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_HPP
#define FOUR_C_PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Teuchos_ParameterList.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
{
  class DEMContactTangentialBase
  {
   public:
    //! constructor
    explicit DEMContactTangentialBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactTangentialBase() = default;

    //! init tangential contact handler
    virtual void init();

    //! setup tangential contact handler
    virtual void setup(const double& k_normal);

    //! set current step size
    virtual void set_current_step_size(const double currentstepsize) final;

    //! calculate tangential contact force
    virtual void tangential_contact_force(double* gap_tangential, bool& stick_tangential,
        const double* normal, const double* v_rel_tangential, const double& m_eff,
        const double& mu_tangential, const double& normalcontactforce,
        double* tangentialcontactforce) const = 0;

    //! evaluate tangential potential energy
    virtual void tangential_potential_energy(
        const double* gap_tangential, double& tangentialpotentialenergy) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! timestep size
    double dt_;
  };

  class DEMContactTangentialLinearSpringDamp : public DEMContactTangentialBase
  {
   public:
    //! constructor
    explicit DEMContactTangentialLinearSpringDamp(const Teuchos::ParameterList& params);

    //! init tangential contact handler
    void init() override;

    //! setup tangential contact handler
    void setup(const double& k_normal) override;

    //! calculate tangential contact force
    void tangential_contact_force(double* gap_tangential, bool& stick_tangential,
        const double* normal, const double* v_rel_tangential, const double& m_eff,
        const double& mu_tangential, const double& normalcontactforce,
        double* tangentialcontactforce) const override;

    //! evaluate tangential potential energy
    void tangential_potential_energy(
        const double* gap_tangential, double& tangentialpotentialenergy) const override;

   private:
    //! coefficient of restitution
    const double e_;

    //! particle Poisson ratio
    const double nue_;

    //! tangential contact stiffness
    double k_tangential_;

    //! tangential contact damping factor
    double d_tangential_fac_;
  };

}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
