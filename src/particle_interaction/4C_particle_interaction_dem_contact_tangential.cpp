/*---------------------------------------------------------------------------*/
/*! \file
\brief tangential contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_dem_contact_tangential.hpp"

#include "4C_inpar_particle.hpp"
#include "4C_particle_interaction_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::DEMContactTangentialBase::DEMContactTangentialBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params), dt_(0.0)
{
  // empty constructor
}

void ParticleInteraction::DEMContactTangentialBase::init()
{
  // nothing to do
}

void ParticleInteraction::DEMContactTangentialBase::setup(const double& k_normal)
{
  // nothing to do
}

void ParticleInteraction::DEMContactTangentialBase::set_current_step_size(
    const double currentstepsize)
{
  dt_ = currentstepsize;
}

ParticleInteraction::DEMContactTangentialLinearSpringDamp::DEMContactTangentialLinearSpringDamp(
    const Teuchos::ParameterList& params)
    : ParticleInteraction::DEMContactTangentialBase(params),
      e_(params_dem_.get<double>("COEFF_RESTITUTION")),
      nue_(params_dem_.get<double>("POISSON_RATIO")),
      k_tangential_(0.0),
      d_tangential_fac_(0.0)
{
  // empty constructor
}

void ParticleInteraction::DEMContactTangentialLinearSpringDamp::init()
{
  // call base class init
  DEMContactTangentialBase::init();

  // safety checks for contact parameters
  if (nue_ <= -1.0 or nue_ > 0.5)
    FOUR_C_THROW("invalid input parameter POISSON_RATIO (expected in range ]-1.0; 0.5])!");

  if (params_dem_.get<double>("FRICT_COEFF_TANG") <= 0.0)
    FOUR_C_THROW("invalid input parameter FRICT_COEFF_TANG for this kind of contact law!");
}

void ParticleInteraction::DEMContactTangentialLinearSpringDamp::setup(const double& k_normal)
{
  // call base class setup
  DEMContactTangentialBase::setup(k_normal);

  // tangential to normal stiffness ratio
  const double kappa = (1.0 - nue_) / (1.0 - 0.5 * nue_);

  // tangential contact stiffness
  k_tangential_ = kappa * k_normal;

  // determine tangential contact damping factor
  if (e_ > 0.0)
  {
    const double lne = std::log(e_);
    d_tangential_fac_ =
        2.0 * std::abs(lne) * std::sqrt(k_normal / (UTILS::Pow<2>(lne) + UTILS::Pow<2>(M_PI)));
  }
  else
    d_tangential_fac_ = 2.0 * std::sqrt(k_normal);
}

void ParticleInteraction::DEMContactTangentialLinearSpringDamp::tangential_contact_force(
    double* gap_tangential, bool& stick_tangential, const double* normal,
    const double* v_rel_tangential, const double& m_eff, const double& mu_tangential,
    const double& normalcontactforce, double* tangentialcontactforce) const
{
  // determine tangential contact damping parameter
  const double d_tangential = d_tangential_fac_ * std::sqrt(m_eff);

  // compute length of tangential gap at time n
  const double old_length = UTILS::VecNormTwo(gap_tangential);

  // compute projection of tangential gap onto current normal at time n+1
  UTILS::VecAddScale(gap_tangential, -UTILS::VecDot(normal, gap_tangential), normal);

  // compute length of tangential gap at time n+1
  const double new_length = UTILS::VecNormTwo(gap_tangential);

  // maintain length of tangential gap equal to before the projection
  if (new_length > 1.0e-14)
    UTILS::VecSetScale(gap_tangential, old_length / new_length, gap_tangential);

  // update of elastic tangential displacement if stick is true
  if (stick_tangential == true) UTILS::VecAddScale(gap_tangential, dt_, v_rel_tangential);

  // compute tangential contact force (assume stick-case)
  UTILS::VecSetScale(tangentialcontactforce, -k_tangential_, gap_tangential);
  UTILS::VecAddScale(tangentialcontactforce, -d_tangential, v_rel_tangential);

  // compute the norm of the tangential contact force
  const double norm_tangentialcontactforce = UTILS::VecNormTwo(tangentialcontactforce);

  // tangential contact force for stick-case
  if (norm_tangentialcontactforce <= (mu_tangential * std::abs(normalcontactforce)))
  {
    stick_tangential = true;

    // tangential contact force already computed
  }
  // tangential contact force for slip-case
  else
  {
    stick_tangential = false;

    // compute tangential contact force
    UTILS::VecSetScale(tangentialcontactforce,
        mu_tangential * std::abs(normalcontactforce) / norm_tangentialcontactforce,
        tangentialcontactforce);

    // compute tangential displacement
    const double inv_k_tangential = 1.0 / k_tangential_;
    UTILS::VecSetScale(gap_tangential, -inv_k_tangential, tangentialcontactforce);
    UTILS::VecAddScale(gap_tangential, -inv_k_tangential * d_tangential, v_rel_tangential);
  }
}

void ParticleInteraction::DEMContactTangentialLinearSpringDamp::tangential_potential_energy(
    const double* gap_tangential, double& tangentialpotentialenergy) const
{
  tangentialpotentialenergy = 0.5 * k_tangential_ * UTILS::VecDot(gap_tangential, gap_tangential);
}

FOUR_C_NAMESPACE_CLOSE
