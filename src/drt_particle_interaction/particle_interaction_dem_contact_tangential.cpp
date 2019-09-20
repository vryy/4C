/*---------------------------------------------------------------------------*/
/*! \file
\brief tangential contact handler for discrete element method (DEM) interactions

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_contact_tangential.H"

#include "particle_interaction_utils.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContactTangentialBase::DEMContactTangentialBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params), dt_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactTangentialBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactTangentialBase::Setup(const double& k_normal)
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactTangentialBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactTangentialBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactTangentialBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp::DEMContactTangentialLinearSpringDamp(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactTangentialBase(params),
      e_(params_dem_.get<double>("COEFF_RESTITUTION")),
      nue_(params_dem_.get<double>("POISSON_RATIO")),
      k_tangential_(0.0),
      d_tangential_fac_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp::Init()
{
  // call base class init
  DEMContactTangentialBase::Init();

  // safety checks for contact parameters
  if (nue_ <= -1.0 or nue_ > 0.5)
    dserror("invalid input parameter POISSON_RATIO (expected in range ]-1.0; 0.5])!");

  if (params_dem_.get<double>("FRICT_COEFF_TANG") <= 0.0)
    dserror("invalid input parameter FRICT_COEFF_TANG for this kind of contact law!");
}

void PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp::Setup(const double& k_normal)
{
  // call base class setup
  DEMContactTangentialBase::Setup(k_normal);

  // tangential to normal stiffness ratio
  const double kappa = (1.0 - nue_) / (1.0 - 0.5 * nue_);

  // tangential contact stiffness
  k_tangential_ = kappa * k_normal;

  // determine tangential contact damping factor
  if (e_ > 0.0)
  {
    const double lne = std::log(e_);
    d_tangential_fac_ =
        2.0 * std::abs(lne) * std::sqrt(k_normal / (UTILS::pow<2>(lne) + UTILS::pow<2>(M_PI)));
  }
  else
    d_tangential_fac_ = 2.0 * std::sqrt(k_normal);
}

void PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp::TangentialContactForce(
    double* gap_tangential, bool& stick_tangential, const double* normal,
    const double* v_rel_tangential, const double& m_eff, const double& mu_tangential,
    const double& normalcontactforce, double* tangentialcontactforce) const
{
  // determine tangential contact damping parameter
  const double d_tangential = d_tangential_fac_ * std::sqrt(m_eff);

  // compute length of tangential gap at time n
  const double old_length = UTILS::vec_norm2(gap_tangential);

  // compute projection of tangential gap onto current normal at time n+1
  UTILS::vec_addscale(gap_tangential, -UTILS::vec_dot(normal, gap_tangential), normal);

  // compute length of tangential gap at time n+1
  const double new_length = UTILS::vec_norm2(gap_tangential);

  // maintain length of tangential gap equal to before the projection
  if (new_length > 1.0e-14)
    UTILS::vec_setscale(gap_tangential, old_length / new_length, gap_tangential);

  // update of elastic tangential displacement if stick is true
  if (stick_tangential == true) UTILS::vec_addscale(gap_tangential, dt_, v_rel_tangential);

  // compute tangential contact force (assume stick-case)
  UTILS::vec_setscale(tangentialcontactforce, -k_tangential_, gap_tangential);
  UTILS::vec_addscale(tangentialcontactforce, -d_tangential, v_rel_tangential);

  // compute the norm of the tangential contact force
  const double norm_tangentialcontactforce = UTILS::vec_norm2(tangentialcontactforce);

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
    UTILS::vec_setscale(tangentialcontactforce,
        mu_tangential * std::abs(normalcontactforce) / norm_tangentialcontactforce,
        tangentialcontactforce);

    // compute tangential displacement
    const double inv_k_tangential = 1.0 / k_tangential_;
    UTILS::vec_setscale(gap_tangential, -inv_k_tangential, tangentialcontactforce);
    UTILS::vec_addscale(gap_tangential, -inv_k_tangential * d_tangential, v_rel_tangential);
  }
}
