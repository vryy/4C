/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_contact_normal.H

\brief normal contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H
#define PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalBase(){};

    //! init normal contact handler
    virtual void Init();

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! write restart of normal contact handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of normal contact handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get normal contact stiffness
    virtual double GetNormalContactStiffness() const = 0;

    //! get time critical contact stiffness
    virtual double GetTimeCriticalStiffness() const = 0;

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! maximum expected particle radius
    const double r_max_;

    //! maximum expected particle velocity
    const double v_max_;

    //! maximum allowed relative penetration (particle-particle)
    const double c_;

    //! normal contact stiffness (particle-particle)
    double k_normal_;
  };

  class DEMContactNormalLinearSpring : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpring(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalLinearSpring(){};

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! get normal contact stiffness
    virtual double GetNormalContactStiffness() const override { return k_normal_; };

    //! get time critical contact stiffness
    virtual double GetTimeCriticalStiffness() const override { return k_normal_; };

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;
  };

  class DEMContactNormalLinearSpringDamp : public DEMContactNormalLinearSpring
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpringDamp(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalLinearSpringDamp(){};

    //! init normal contact handler
    virtual void Init();

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;

   private:
    //! coefficient of restitution (particle-particle)
    const double e_;

    //! parameter for regularization of damping contact force
    const double damp_reg_fac_;

    //! tension cutoff
    const bool tension_cutoff_;

    //! normal contact damping factor (particle-particle)
    double d_normal_fac_;
  };

  class DEMContactNormalNonlinearBase : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalNonlinearBase(){};

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! get normal contact stiffness
    virtual double GetNormalContactStiffness() const override { return k_normal_; };

    //! get time critical contact stiffness
    virtual double GetTimeCriticalStiffness() const override { return k_tcrit_; };

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const = 0;

   private:
    //! time step critical contact stiffness (particle-particle)
    double k_tcrit_;
  };

  class DEMContactNormalHertz : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalHertz(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalHertz(){};

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;
  };

  class DEMContactNormalNonlinearDampBase : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearDampBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalNonlinearDampBase(){};

    //! init normal contact handler
    virtual void Init();

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const = 0;

   protected:
    //! normal contact damping parameter (particle-particle)
    const double d_normal_;

    //! tension cutoff
    const bool tension_cutoff_;
  };

  class DEMContactNormalLeeHerrmann : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLeeHerrmann(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalLeeHerrmann(){};

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;
  };

  class DEMContactNormalKuwabaraKono : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalKuwabaraKono(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalKuwabaraKono(){};

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;
  };

  class DEMContactNormalTsuji : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalTsuji(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalTsuji(){};

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double* radius_i,
        const double* radius_j, const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
