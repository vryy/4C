/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_contact_tangential.H

\brief tangential contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_H
#define PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactTangentialBase
  {
   public:
    //! constructor
    explicit DEMContactTangentialBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactTangentialBase(){};

    //! init tangential contact handler
    virtual void Init();

    //! setup tangential contact handler
    virtual void Setup(const double& k_normal);

    //! write restart of tangential contact handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of tangential contact handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize) final;

    //! calculate tangential contact force
    virtual void TangentialContactForce(double* gap_tangential, bool& stick_tangential,
        const double* normal, const double* v_rel_tangential, const double& m_eff,
        const double& normalcontactforce, double* tangentialcontactforce) const = 0;

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

    //! virtual destructor
    virtual ~DEMContactTangentialLinearSpringDamp(){};

    //! init tangential contact handler
    virtual void Init();

    //! setup tangential contact handler
    virtual void Setup(const double& k_normal);

    //! calculate tangential contact force
    virtual void TangentialContactForce(double* gap_tangential, bool& stick_tangential,
        const double* normal, const double* v_rel_tangential, const double& m_eff,
        const double& normalcontactforce, double* tangentialcontactforce) const override;

   private:
    //! coefficient of restitution (particle-particle)
    const double e_;

    //! particle Poisson ratio
    const double nue_;

    //! tangential contact friction coefficient (particle-particle)
    const double mu_tangential_;

    //! tangential contact stiffness (particle-particle)
    double k_tangential_;

    //! tangential contact damping factor (particle-particle)
    double d_tangential_fac_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
