/*--------------------------------------------------------------------------*/
/*!
\file particle_timint_strategy.cpp

\brief time integration strategies for particle problems

\level 1

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "particle_timint_strategy.H"

#include "particle_algorithm.H"
#include "particle_utils.H"

#include "particle_ellipsoid_node.H"
#include "particle_radius_node.H"
#include "particle_timint.H"

#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/particle_mat_ellipsoids.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::TimIntStrategyBase::TimIntStrategyBase(TimInt* const timint  //!< time integrator
    )
    : timint_(timint)
{
  return;
}


/*----------------------------------------------------------------------*
 | constructor                                               fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::TimIntStrategySpheres::TimIntStrategySpheres(TimInt* const timint  //!< time integrator
    )
    :  // call base class constructor
      TimIntStrategyBase(timint)
{
  // initialize radius vector
  timint->radius_ =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, timint->NodeRowMapView(), true));

  return;
};


/*----------------------------------------------------------------------*
 | compute angular acceleration vector                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::ComputeAngularAcceleration(
    Epetra_Vector& angacc,          //!< global angular acceleration vector
    const Epetra_Vector& m_contact  //!< global moment vector
    ) const
{
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    const double invinertia = 1. / (*timint_->inertia_)[i];
    for (int dim = 0; dim < 3; ++dim)
    {
      angacc[i * 3 + dim] = invinertia * m_contact[i * 3 + dim];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute inertia vector                                    fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::ComputeInertia() const
{
  // safety checks
  if ((*timint_->radius_)(0) == Teuchos::null or timint_->mass_ == Teuchos::null)
    dserror("Invalid vector!");

  // initialize inertia vector
  timint_->inertia_ = Teuchos::rcp(new Epetra_Vector(timint_->mass_->Map(), true));

  // compute inertia for every spherical particle: I = 2/5 * m * r^2
  for (int inode = 0; inode < timint_->mass_->MyLength(); ++inode)
  {
    const double radius = (*(*timint_->radius_)(0))[inode];
    (*timint_->inertia_)[inode] = .4 * (*timint_->mass_)[inode] * radius * radius;
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute kinetic energy                                    fang 10/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::TimIntStrategySpheres::ComputeKineticEnergy() const
{
  double kineticenergy(0.);

  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    double kineticenergy_trans = 0.;
    double kineticenergy_rot = 0.;

    for (int dim = 0; dim < 3; ++dim)
    {
      // translational kinetic energy
      kineticenergy_trans += (*timint_->veln_)[i * 3 + dim] * (*timint_->veln_)[i * 3 + dim];

      // rotational kinetic energy
      kineticenergy_rot += (*timint_->angVeln_)[i * 3 + dim] * (*timint_->angVeln_)[i * 3 + dim];
    }

    kineticenergy += .5 * ((*timint_->mass_)[i] * kineticenergy_trans +
                              (*timint_->inertia_)[i] * kineticenergy_rot);
  }

  return kineticenergy;
}


/*----------------------------------------------------------------------*
 | compute mass vector                                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::ComputeMass() const
{
  // loop over all particles
  for (int i = 0; i < timint_->discret_->NumMyRowNodes(); ++i)
    // particle mass: m = rho * 4/3 * PI *r^3
    (*timint_->mass_)[i] = timint_->particle_algorithm_->ParticleMat()[0]->initDensity_ *
                           PARTICLE::Utils::Radius2Volume((*(*timint_->radius_)(0))[i]);

  return;
}


/*----------------------------------------------------------------------*
 | output particle orientation                               fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::OutputOrientation() const
{
  if (timint_->writeorientation_) timint_->output_->WriteVector("orientation", timint_->orient_);

  return;
}


/*----------------------------------------------------------------------*
 | predict or correct angular velocity vector                fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::PredictOrCorrectAngularVelocity(
    const double dt,                 //!< time step size
    const Epetra_Vector& m_contact,  //!< global moment vector
    const Epetra_Vector& orient      //!< particle orientation vector at time t_n
    ) const
{
  // output warning
  dserror(
      "Prediction and correction of angular velocity not yet tested for spherical particles! "
      "Delete this error if you want, but make sure that the following code does the right job...");

  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    const double invinertia = 1. / (*timint_->inertia_)[i];
    for (int dim = 0; dim < 3; ++dim)
      (*timint_->angVeln_)[i * 3 + dim] += invinertia * dt * (m_contact)[i * 3 + dim];
  }

  return;
}


/*----------------------------------------------------------------------*
 | update particle orientation vector                        fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::RotateOrientVector(const double dt  //!< time step size
    ) const
{
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    double angVel[3];
    double r[3];

    for (int dim = 0; dim < 3; ++dim)
    {
      angVel[dim] = (*timint_->angVeln_)[i * 3 + dim];
      r[dim] = (*timint_->orient_)[i * 3 + dim];
    }

    // visualization only valid for 2D when rotating around y-axis

    // simplified/linearized orient vector - just for visualization
    // delta_r = \Delta t * (ang_vel x r)
    (*timint_->orient_)[i * 3] += dt * (angVel[1] * r[2] - angVel[2] * r[1]);
    (*timint_->orient_)[i * 3 + 1] += dt * (angVel[2] * r[0] - angVel[0] * r[2]);
    (*timint_->orient_)[i * 3 + 2] += dt * (angVel[0] * r[1] - angVel[1] * r[0]);
    //--------------------------------------------------------------

    // more exactly------------------------------------------------
    //    double d[3];
    //    double norm_ang_vel=0.0;
    //    //norm ang_vel
    //    for(int dim=0; dim<3; ++dim)
    //    {
    //      norm_ang_vel += ang_vel[dim]*ang_vel[dim];
    //    }
    //    norm_ang_vel = sqrt(norm_ang_vel);
    //
    //    if(norm_ang_vel > 1E-10)
    //    {
    //      double scalar = 0.0;
    //      double invnorm_ang_vel = 1.0 / norm_ang_vel;
    //      for(int dim=0; dim<3; ++dim)
    //      {
    //        d[dim] = invnorm_ang_vel * ang_vel[dim];
    //        scalar += d[dim] * r[dim];
    //      }
    //
    //      for(int dim=0; dim<3; ++dim)
    //      {
    //        (*orient_)[i*3+dim] += (1-cos(norm_ang_vel*dt)) * scalar * d[dim] -
    //        (1-cos(norm_ang_vel*dt)) * r[dim]
    //              + sin(norm_ang_vel*dt) * ( d[(dim+1)%3] * r[(dim+2)%3] - d[(dim+2)%3] *
    //              r[(dim+1)%3] );
    //      }
    //    }
    //--------------------------------------------------------------
  }

  return;
}


/*----------------------------------------------------------------------*
 | set initial particle orientation                          fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::SetInitialOrientation() const
{
  dserror("Function not yet implemented for spherical particles!");

  return;
}


/*----------------------------------------------------------------------*
 | set initial particle radii                                fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::SetInitialRadii() const
{
  // extract initial particle radius if available
  std::vector<const MAT::PAR::ParticleMat*> particlemat =
      timint_->particle_algorithm_->ParticleMat();
  // safety check
  if (particlemat.size() == 2 and particlemat[0]->initRadius_ != particlemat[1]->initRadius_)
    dserror("Both particle materials need to be defined with equal initRadius_!");

  const double initRadius = particlemat.size() > 0 ? particlemat[0]->initRadius_ : 0.0;

  // set initial particle radii to initial or specified values
  for (int i = 0; i < timint_->discret_->NumMyRowNodes(); ++i)
  {
    const ParticleRadiusNode* const particle =
        dynamic_cast<const ParticleRadiusNode* const>(timint_->discret_->lRowNode(i));
    (*(*timint_->radius_)(0))[i] = particle == NULL ? initRadius : particle->Radius();
  }

  return;
}


/*----------------------------------------------------------------------*
 | update map of inertia vector                              fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::UpdateInertiaVectorMap() const
{
  timint_->UpdateStateVectorMap(timint_->inertia_, true);

  return;
}


/*----------------------------------------------------------------------*
 | update map of radius vector                               fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategySpheres::UpdateRadiusVectorMap() const
{
  timint_->UpdateStateVectorMap(timint_->radius_, true);

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                               fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::TimIntStrategyEllipsoids::TimIntStrategyEllipsoids(
    TimInt* const timint  //!< time integrator
    )
    :  // call base class constructor
      TimIntStrategyBase(timint)
{
  // safety checks
  if (DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(
          DRT::Problem::Instance()->ParticleParams(), "DYNAMICTYP") !=
          INPAR::PARTICLE::dyna_centrdiff and
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(
          DRT::Problem::Instance()->ParticleParams(), "DYNAMICTYP") != INPAR::PARTICLE::dyna_rk2)
    dserror("Only central differencing and RK2 available for ellipsoids!");

  // set flag in time integrator
  timint_->writeorientation_ = true;

  // initialize radius vector to store semi-axes of ellipsoidal particles
  timint->radius_ =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, timint->DofRowMapView(), true));

  return;
};


/*----------------------------------------------------------------------*
 | compute angular acceleration vector                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::ComputeAngularAcceleration(
    Epetra_Vector& angacc,          //!< global angular acceleration vector
    const Epetra_Vector& m_contact  //!< global moment vector
    ) const
{
  // initialize vectors for calculation
  static LINALG::Matrix<3, 1> inertia, orient, angVeln, moment, angAccn, temp;

  // initialize rotation matrix
  static LINALG::Matrix<3, 3> R;

  // loop over all particles
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    // extract variables associated with current particle
    for (int dim = 0; dim < 3; ++dim)
    {
      inertia(dim) = (*timint_->inertia_)[i * 3 + dim];
      orient(dim) = (*timint_->orient_)[i * 3 + dim];
      angVeln(dim) = (*timint_->angVeln_)[i * 3 + dim];
      moment(dim) = (m_contact)[i * 3 + dim];
    }

    // compute rotation matrix
    LARGEROTATIONS::angletotriad(orient, R);

    // transform angular velocity and moment to particle coordinate system
    temp.Multiply(R, angVeln);
    angVeln = temp;
    temp.Multiply(R, moment);
    moment = temp;

    // compute angular acceleration expressed in particle coordinate system:
    // angAccn = inertia^-1 * (moment - angVeln x (inertia * angVeln))
    temp.EMultiply(inertia, angVeln);
    angAccn.CrossProduct(angVeln, temp);
    moment -= angAccn;
    moment.EDivide(inertia);

    // transform angular acceleration back to global coordinate system
    angAccn.MultiplyTN(R, moment);

    // insert angular acceleration of current particle into global vector
    for (int dim = 0; dim < 3; ++dim) angacc[i * 3 + dim] = angAccn(dim);
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute inertia vector                                    fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::ComputeInertia() const
{
  // safety check
  if (timint_->mass_ == Teuchos::null) dserror("Invalid mass vector!");

  // initialize inertia vector
  timint_->inertia_ = Teuchos::rcp(new Epetra_Vector(*timint_->discret_->DofRowMap(), true));

  // compute inertia for every particle
  for (int inode = 0; inode < timint_->mass_->MyLength(); ++inode)
  {
    // precompute quantities
    const double prefactor = .2 * (*timint_->mass_)[inode];
    const double semiaxis_0 = (*(*timint_->radius_)(0))[inode * 3];
    const double semiaxis_1 = (*(*timint_->radius_)(0))[inode * 3 + 1];
    const double semiaxis_2 = (*(*timint_->radius_)(0))[inode * 3 + 2];
    const double semiaxis_0_sq = semiaxis_0 * semiaxis_0;
    const double semiaxis_1_sq = semiaxis_1 * semiaxis_1;
    const double semiaxis_2_sq = semiaxis_2 * semiaxis_2;
    (*timint_->inertia_)[inode * 3] = prefactor * (semiaxis_1_sq + semiaxis_2_sq);
    (*timint_->inertia_)[inode * 3 + 1] = prefactor * (semiaxis_0_sq + semiaxis_2_sq);
    (*timint_->inertia_)[inode * 3 + 2] = prefactor * (semiaxis_0_sq + semiaxis_1_sq);
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute kinetic energy                                    fang 10/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::TimIntStrategyEllipsoids::ComputeKineticEnergy() const
{
  // initialize result variable
  double kineticenergy(0.);

  // initialize vectors for calculation
  static LINALG::Matrix<3, 1> inertia, orient, angVeln, temp;

  // initialize rotation matrix
  static LINALG::Matrix<3, 3> R;

  // loop over all particles
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    // translational kinetic energy
    double kineticenergy_trans = 0.;
    for (int dim = 0; dim < 3; ++dim)
    {
      const double veln = (*timint_->veln_)[i * 3 + dim];
      kineticenergy_trans += veln * veln;
    }
    kineticenergy_trans *= .5 * (*timint_->mass_)[i];

    // rotational kinetic energy
    for (int dim = 0; dim < 3; ++dim)
    {
      inertia(dim) = (*timint_->inertia_)[i * 3 + dim];
      orient(dim) = (*timint_->orient_)[i * 3 + dim];
      angVeln(dim) = (*timint_->angVeln_)[i * 3 + dim];
    }
    LARGEROTATIONS::angletotriad(orient, R);
    temp.Multiply(R, angVeln);
    angVeln = temp;
    inertia.EMultiply(angVeln);
    const double kineticenergy_rot = .5 * inertia.Dot(angVeln);

    // sum up translational and rotational contributions
    kineticenergy += kineticenergy_trans + kineticenergy_rot;
  }

  return kineticenergy;
}


/*----------------------------------------------------------------------*
 | compute mass vector                                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::ComputeMass() const
{
  // loop over all particles
  for (int i = 0; i < timint_->discret_->NumMyRowNodes(); ++i)
    // particle mass: m = rho * 4/3 * PI * r1 * r2 * r3
    (*timint_->mass_)[i] = timint_->particle_algorithm_->ParticleMat()[0]->initDensity_ * 4. / 3. *
                           M_PI * (*(*timint_->radius_)(0))[i * 3] *
                           (*(*timint_->radius_)(0))[i * 3 + 1] *
                           (*(*timint_->radius_)(0))[i * 3 + 2];

  return;
}


/*----------------------------------------------------------------------*
 | output particle orientation                               fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::OutputOrientation() const
{
  // output standard orientation vector
  timint_->output_->WriteVector("orientation", timint_->orient_);

  // initialize vector for storage of orientation tensors associated with single particles
  const Teuchos::RCP<Epetra_MultiVector> orient_tensor =
      Teuchos::rcp(new Epetra_MultiVector(*timint_->discret_->NodeRowMap(), 9));

  // initialize vector and matrices for calculation of orientation tensors
  static LINALG::Matrix<3, 1> orient;
  static LINALG::Matrix<3, 3> E, E_rot, R, temp;

  // loop over all particles
  for (int inode = 0; inode < timint_->discret_->NodeRowMap()->NumMyElements(); ++inode)
  {
    // construct plain ellipsoid matrix without considering translation and rotation
    E.Clear();
    for (unsigned dim = 0; dim < 3; ++dim) E(dim, dim) = (*(*timint_->radius_)(0))[inode * 3 + dim];

    // extract orientation of current particle
    for (int dim = 0; dim < 3; ++dim) orient(dim) = (*timint_->orient_)[inode * 3 + dim];

    // compute rotation matrix for current particle
    LARGEROTATIONS::angletotriad(orient, R);

    // apply rotation matrix to ellipsoid matrix
    temp.Multiply(R, E);
    E_rot.MultiplyNT(temp, R);

    // write final orientation tensor into vector
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) orient_tensor->ReplaceMyValue(inode, j * 3 + i, E_rot(i, j));
  }

  // output vector containing orientation tensors associated with single particles
  timint_->output_->WriteVector("orientation_tensor", orient_tensor);

  return;
}


/*----------------------------------------------------------------------*
 | predict or correct angular velocity vector                fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::PredictOrCorrectAngularVelocity(
    const double dt,                 //!< time step size
    const Epetra_Vector& m_contact,  //!< global moment vector
    const Epetra_Vector& orient      //!< particle orientation vector at time t_n
    ) const
{
  // initialize vectors and rotation matrix
  static LINALG::Matrix<3, 1> angVeln, inertia, moment, orientold, orientn, temp;
  static LINALG::Matrix<3, 3> R;

  // loop over all particles
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    // extract variables associated with current particle
    for (unsigned dim = 0; dim < 3; ++dim)
    {
      inertia(dim) = (*timint_->inertia_)[i * 3 + dim];
      orientold(dim) = orient[i * 3 + dim];
      orientn(dim) = (*timint_->orient_)[i * 3 + dim];
      angVeln(dim) = (*timint_->angVeln_)[i * 3 + dim];
      moment(dim) = (m_contact)[i * 3 + dim];
    }

    // compute rotation matrix associated with time t_n
    LARGEROTATIONS::angletotriad(orientold, R);

    // transform angular velocity at time t_n to particle coordinate system
    temp.Multiply(R, angVeln);
    angVeln = temp;

    // multiply transformed angular velocity at time t_n by inertia vector to obtain angular
    // momentum at time t_n
    angVeln.EMultiply(inertia);

    // transform angular momentum at time t_n back to global coordinate system
    temp.MultiplyTN(R, angVeln);

    // add change in angular momentum to obtain angular momentum at new time
    temp.Update(dt, moment, 1.);

    // compute rotation matrix associated with new time
    LARGEROTATIONS::angletotriad(orientn, R);

    // transform angular momentum at new time to particle coordinate system
    angVeln.Multiply(R, temp);

    // divide transformed angular momentum at new time by inertia vector to obtain angular velocity
    // at new time
    angVeln.EDivide(inertia);

    // transform angular velocity at new time back to global coordinate system
    temp.MultiplyTN(R, angVeln);

    // insert angular velocity at new time into global vector
    for (unsigned dim = 0; dim < 3; ++dim) (*timint_->angVeln_)[i * 3 + dim] = temp(dim);
  }

  return;
}


/*----------------------------------------------------------------------*
 | update particle orientation vector                        fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::RotateOrientVector(const double dt  //!< time step size
    ) const
{
  // initialize vectors
  static LINALG::Matrix<3, 1> angVeln, orient;
  static LINALG::Matrix<4, 1> q, qinc, qn;

  // loop over all particles
  for (int i = 0; i < timint_->discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    // extract variables associated with current particle
    for (unsigned dim = 0; dim < 3; ++dim)
    {
      orient(dim) = (*timint_->orient_)[i * 3 + dim];
      angVeln(dim) = (*timint_->angVeln_)[i * 3 + dim];
    }

    // compute angular velocity
    const double angVel = angVeln.Norm2();

    // update orientation vector in case of non-vanishing angular velocity
    if (angVel > 1.e-16)
    {
      // express old orientation as quaternion
      LARGEROTATIONS::angletoquaternion(orient, q);

      // compute inverse of angular velocity
      const double angVelinv = 1. / angVel;

      // compute half of the rotation angle increment and corresponding sine
      const double phi = .5 * dt * angVel;
      const double sinphi = sin(phi);

      // express rotation angle increment as quaternion
      qinc(0) = sinphi * angVeln(0) * angVelinv;
      qinc(1) = sinphi * angVeln(1) * angVelinv;
      qinc(2) = sinphi * angVeln(2) * angVelinv;
      qinc(3) = cos(phi);

      // add rotation angle increment via quaternion product
      LARGEROTATIONS::quaternionproduct(q, qinc, qn);

      // express new orientation as rotation angle
      LARGEROTATIONS::quaterniontoangle(qn, orient);

      // insert new orientation into global vector
      for (int dim = 0; dim < 3; ++dim) (*timint_->orient_)[i * 3 + dim] = orient(dim);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | set initial particle orientation                          fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::SetInitialOrientation() const
{
  // safety check
  if (timint_->orient_ == Teuchos::null)
    dserror("Particle orientation vector has not yet been initialized!");

  // loop over all particles
  for (int i = 0; i < timint_->discret_->NumMyRowNodes(); ++i)
  {
    // extract current particle
    const ParticleEllipsoidNode* const particle =
        dynamic_cast<const ParticleEllipsoidNode* const>(timint_->discret_->lRowNode(i));
    if (particle == NULL) dserror("Couldn't extract ellipsoidal particle!");

    // set initial orientation of current particle
    for (unsigned dim = 0; dim < 3; ++dim)
      (*timint_->orient_)[i * 3 + dim] = particle->Orientation()(dim);
  }

  return;
}


/*----------------------------------------------------------------------*
 | set initial particle radii                                fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::SetInitialRadii() const
{
  // extract material parameters for ellipsoidal particles
  const int matid =
      DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids);
  if (matid < 0) dserror("Invalid material ID!");
  const MAT::PAR::ParticleMatEllipsoids* const particlematellipsoids =
      dynamic_cast<const MAT::PAR::ParticleMatEllipsoids* const>(
          DRT::Problem::Instance()->Materials()->ParameterById(matid));
  if (particlematellipsoids == NULL)
    dserror("Couldn't extract material parameters for ellipsoidal particles!");

  // initialize pointer to semi-axes of ellipsoidal particles
  const LINALG::Matrix<3, 1>* semiaxes(NULL);

  // set initial particle semi-axes to default or individually specified values
  for (int i = 0; i < timint_->discret_->NumMyRowNodes(); ++i)
  {
    // extract current particle
    const ParticleEllipsoidNode* const particle =
        dynamic_cast<const ParticleEllipsoidNode* const>(timint_->discret_->lRowNode(i));
    if (particle == NULL) dserror("Couldn't extract ellipsoidal particle!");

    // first case: individually specified semi-axes
    if (particle->SemiAxes()(0) > 0. and particle->SemiAxes()(1) > 0. and
        particle->SemiAxes()(2) > 0.)
      semiaxes = &particle->SemiAxes();

    // second case: default semi-axes
    else if (abs(particle->SemiAxes()(0)) < 1.e-16 and abs(particle->SemiAxes()(1)) < 1.e-16 and
             abs(particle->SemiAxes()(2)) < 1.e-16 and particlematellipsoids->SemiAxes()(0) > 0. and
             particlematellipsoids->SemiAxes()(1) > 0. and
             particlematellipsoids->SemiAxes()(2) > 0.)
      semiaxes = &particlematellipsoids->SemiAxes();

    // third case: invalid semi-axes
    else
      dserror("Must have positive semi-axes!");

    // set semi-axes
    for (unsigned dim = 0; dim < 3; ++dim)
      (*(*timint_->radius_)(0))[i * 3 + dim] = (*semiaxes)(dim);
  }

  return;
}


/*----------------------------------------------------------------------*
 | update map of inertia vector                              fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::UpdateInertiaVectorMap() const
{
  timint_->UpdateStateVectorMap(timint_->inertia_);

  return;
}


/*----------------------------------------------------------------------*
 | update map of radius vector                               fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntStrategyEllipsoids::UpdateRadiusVectorMap() const
{
  timint_->UpdateStateVectorMap(timint_->radius_);

  return;
}
