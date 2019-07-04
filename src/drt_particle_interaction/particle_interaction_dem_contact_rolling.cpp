/*---------------------------------------------------------------------------*/
/*!
\brief rolling contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_contact_rolling.H"

#include "particle_interaction_utils.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContactRollingBase::DEMContactRollingBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params),
      dt_(0.0),
      e_(params_dem_.get<double>("COEFF_RESTITUTION")),
      nue_(params_dem_.get<double>("POISSON_RATIO")),
      mu_rolling_(params_dem_.get<double>("FRICT_COEFF_ROLL")),
      d_rolling_fac_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init rolling contact handler                               sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingBase::Init()
{
  // safety checks for particle-particle contact parameters
  if (nue_ <= -1.0 or nue_ > 0.5)
    dserror("invalid input parameter POISSON_RATIO (expected in range ]-1.0; 0.5])!");

  if (mu_rolling_ <= 0.0)
    dserror("invalid input parameter FRICT_COEFF_ROLL for this kind of contact law!");
}

/*---------------------------------------------------------------------------*
 | setup rolling contact handler                              sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingBase::Setup(const double& k_normal)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of rolling contact handler                   sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of rolling contact handler                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContactRollingViscous::DEMContactRollingViscous(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactRollingBase(params),
      young_(params_dem_.get<double>("YOUNG_MODULUS")),
      v_max_(params_dem_.get<double>("MAX_VELOCITY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | setup rolling contact handler                              sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingViscous::Setup(const double& k_normal)
{
  // call base class setup
  DEMContactRollingBase::Setup(k_normal);

  // determine rolling contact damping factor
  const double fac = young_ / (1.0 - UTILS::pow<2>(nue_));
  const double c_1 = 1.15344;
  d_rolling_fac_ = mu_rolling_ * (1.0 - e_) / (c_1 * std::pow(fac, 0.4) * std::pow(v_max_, 0.2));
}

/*---------------------------------------------------------------------------*
 | calculate effective radius                                 sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingViscous::EffectiveRadiusParticle(
    const double* radius_i, const double* radius_j, const double& gap, double& r_eff) const
{
  if (radius_j)
    r_eff = (radius_i[0] * radius_j[0]) / (radius_i[0] + radius_j[0]);
  else
    r_eff = radius_i[0];
}

/*---------------------------------------------------------------------------*
 | calculate relative rolling velocity                        sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingViscous::RelativeRollingVelocity(const double& r_eff,
    const double* normal, const double* angvel_i, const double* angvel_j,
    double* v_rel_rolling) const
{
  UTILS::vec_setcross(v_rel_rolling, angvel_i, normal);
  if (angvel_j) UTILS::vec_addcross(v_rel_rolling, normal, angvel_j);
}

/*---------------------------------------------------------------------------*
 | calculate rolling contact moment                           sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingViscous::RollingContactMoment(double* gap_rolling,
    bool& stick_rolling, const double* normal, const double* v_rel_rolling, const double& m_eff,
    const double& r_eff, const double& normalcontactforce, double* rollingcontactmoment) const
{
  // determine rolling contact damping parameter
  const double d_rolling = d_rolling_fac_ * std::pow(0.5 * r_eff, -0.2);

  // compute rolling contact force
  double rollingcontactforce[3];
  UTILS::vec_setscale(rollingcontactforce, -(d_rolling * normalcontactforce), v_rel_rolling);

  // compute rolling contact moment
  UTILS::vec_setcross(rollingcontactmoment, rollingcontactforce, normal);
  UTILS::vec_scale(rollingcontactmoment, r_eff);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContactRollingCoulomb::DEMContactRollingCoulomb(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactRollingBase(params), k_rolling_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | setup rolling contact handler                              sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingCoulomb::Setup(const double& k_normal)
{
  // call base class setup
  DEMContactRollingBase::Setup(k_normal);

  // rolling to normal stiffness ratio
  const double kappa = (1.0 - nue_) / (1.0 - 0.5 * nue_);

  // rolling contact stiffness
  k_rolling_ = kappa * k_normal;

  // determine rolling contact damping factor
  if (e_ > 0.0)
  {
    const double lne = std::log(e_);
    d_rolling_fac_ =
        2.0 * std::abs(lne) * std::sqrt(k_normal / (UTILS::pow<2>(lne) + UTILS::pow<2>(M_PI)));
  }
  else
    d_rolling_fac_ = 2.0 * std::sqrt(k_normal);
}

/*---------------------------------------------------------------------------*
 | calculate effective radius for particle-particle contact   sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingCoulomb::EffectiveRadiusParticle(
    const double* radius_i, const double* radius_j, const double& gap, double& r_eff) const
{
  if (radius_j)
    r_eff =
        ((radius_i[0] + 0.5 * gap) * (radius_j[0] + 0.5 * gap)) / (radius_i[0] + radius_j[0] + gap);
  else
    r_eff = radius_i[0] + gap;
}

/*---------------------------------------------------------------------------*
 | calculate relative rolling velocity                        sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingCoulomb::RelativeRollingVelocity(const double& r_eff,
    const double* normal, const double* angvel_i, const double* angvel_j,
    double* v_rel_rolling) const
{
  UTILS::vec_setcross(v_rel_rolling, normal, angvel_i);
  if (angvel_j) UTILS::vec_addcross(v_rel_rolling, angvel_j, normal);

  UTILS::vec_scale(v_rel_rolling, r_eff);
}

/*---------------------------------------------------------------------------*
 | calculate rolling contact moment                           sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContactRollingCoulomb::RollingContactMoment(double* gap_rolling,
    bool& stick_rolling, const double* normal, const double* v_rel_rolling, const double& m_eff,
    const double& r_eff, const double& normalcontactforce, double* rollingcontactmoment) const
{
  // determine rolling contact damping parameter
  const double d_rolling = d_rolling_fac_ * std::sqrt(m_eff);

  // compute length of rolling gap at time n
  const double old_length = UTILS::vec_norm2(gap_rolling);

  // compute projection of rolling gap onto current normal at time n+1
  UTILS::vec_addscale(gap_rolling, -UTILS::vec_dot(normal, gap_rolling), normal);

  // compute length of rolling gap at time n+1
  const double new_length = UTILS::vec_norm2(gap_rolling);

  // maintain length of rolling gap equal to before the projection
  if (new_length > 1.0e-14) UTILS::vec_setscale(gap_rolling, old_length / new_length, gap_rolling);

  // update of elastic rolling displacement if stick is true
  if (stick_rolling == true) UTILS::vec_addscale(gap_rolling, dt_, v_rel_rolling);

  // compute rolling contact force (assume stick-case)
  double rollingcontactforce[3];
  UTILS::vec_setscale(rollingcontactforce, -k_rolling_, gap_rolling);
  UTILS::vec_addscale(rollingcontactforce, -d_rolling, v_rel_rolling);

  // compute the norm of the rolling contact force
  const double norm_rollingcontactforce = UTILS::vec_norm2(rollingcontactforce);

  // rolling contact force for stick-case
  if (norm_rollingcontactforce <= (mu_rolling_ * std::abs(normalcontactforce)))
  {
    stick_rolling = true;

    // rolling contact force already computed
  }
  // rolling contact force for slip-case
  else
  {
    stick_rolling = false;

    // compute rolling contact force
    UTILS::vec_setscale(rollingcontactforce,
        mu_rolling_ * std::abs(normalcontactforce) / norm_rollingcontactforce, rollingcontactforce);

    // compute rolling displacement
    const double inv_k_rolling = 1.0 / k_rolling_;
    UTILS::vec_setscale(gap_rolling, -inv_k_rolling, rollingcontactforce);
    UTILS::vec_addscale(gap_rolling, -inv_k_rolling * d_rolling, v_rel_rolling);
  }

  // compute rolling contact moment
  UTILS::vec_setcross(rollingcontactmoment, rollingcontactforce, normal);
  UTILS::vec_scale(rollingcontactmoment, r_eff);
}
