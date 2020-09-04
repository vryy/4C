/*---------------------------------------------------------------------------*/
/*! \file
\brief normal contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H
#define PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalBase() = default;

    //! init normal contact handler
    virtual void Init();

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! get normal contact stiffness
    virtual double GetNormalContactStiffness() const = 0;

    //! get time critical contact stiffness
    virtual double GetTimeCriticalStiffness() const = 0;

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const = 0;

    //! evaluate normal potential energy
    virtual void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const = 0;

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
  };

  class DEMContactNormalLinearSpring : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpring(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! get normal contact stiffness
    double GetNormalContactStiffness() const override { return k_normal_; };

    //! get time critical contact stiffness
    double GetTimeCriticalStiffness() const override { return k_normal_; };

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;
  };

  class DEMContactNormalLinearSpringDamp : public DEMContactNormalLinearSpring
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpringDamp(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

   private:
    //! coefficient of restitution
    const double e_;

    //! parameter for regularization of damping contact force
    const double damp_reg_fac_;

    //! tension cutoff
    const bool tension_cutoff_;

    //! normal contact damping factor
    double d_normal_fac_;
  };

  class DEMContactNormalNonlinearBase : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearBase(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! get normal contact stiffness
    double GetNormalContactStiffness() const override { return k_normal_; };

    //! get time critical contact stiffness
    double GetTimeCriticalStiffness() const override { return k_tcrit_; };

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;

   private:
    //! time step critical contact stiffness
    double k_tcrit_;
  };

  class DEMContactNormalHertz : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalHertz(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalNonlinearDampBase : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearDampBase(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

   protected:
    //! normal contact damping parameter
    const double d_normal_;

    //! tension cutoff
    const bool tension_cutoff_;
  };

  class DEMContactNormalLeeHerrmann : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLeeHerrmann(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalKuwabaraKono : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalKuwabaraKono(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalTsuji : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalTsuji(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double* radius_i, const double* radius_j,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
