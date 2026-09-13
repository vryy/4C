// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_HPP
#define FOUR_C_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace Particle
{
  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalBase() = default;

    //! setup normal contact handler
    virtual void setup(const double& dens_max);

    //! get normal contact stiffness
    virtual double get_normal_contact_stiffness() const final { return k_normal_; };

    //! get critical normal contact stiffness
    virtual double get_critical_normal_contact_stiffness() const final { return k_normal_crit_; };

    //! evaluate normal contact force
    virtual void normal_contact_force(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const = 0;

    //! evaluate normal potential energy
    virtual void normal_potential_energy(
        const double& gap, double& normalpotentialenergy) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! maximum expected particle radius
    const double r_max_;

    //! maximum expected particle velocity
    const double v_max_;

    //! maximum allowed relative penetration
    const double c_;

    //! normal contact stiffness
    double k_normal_;

    //! critical normal contact stiffness
    double k_normal_crit_;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalLinearSpring : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpring(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void setup(const double& dens_max) override;

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

    //! evaluate normal potential energy
    void normal_potential_energy(const double& gap, double& normalpotentialenergy) const override;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalLinearSpringDamp
      : public DEMContactNormalLinearSpring
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpringDamp(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void setup(const double& dens_max) override;

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

   private:
    //! coefficient of restitution
    const double e_;

    //! parameter for regularization of damping contact force
    const double damp_reg_fac_;

    //! normal contact damping factor
    double d_normal_fac_;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalNonlinearBase : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearBase(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void setup(const double& dens_max) override;

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

    //! evaluate normal potential energy
    void normal_potential_energy(const double& gap, double& normalpotentialenergy) const override;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalHertz : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalHertz(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalNonlinearDampBase
      : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearDampBase(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

   protected:
    //! normal contact damping parameter
    const double d_normal_;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalLeeHerrmann
      : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLeeHerrmann(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalKuwabaraKono
      : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalKuwabaraKono(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class FOUR_C_API(FOUR_C_CORE) DEMContactNormalTsuji : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalTsuji(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void normal_contact_force(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

}  // namespace Particle

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
