/*--------------------------------------------------------------------------*/
/*!
\file particle_timint_strategy.H

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
#ifndef PARTICLE_TIMINT_STRATEGY_H
#define PARTICLE_TIMINT_STRATEGY_H

#include "../linalg/linalg_fixedsizematrix.H"

// forward declaration
class Epetra_Vector;

namespace PARTICLE
{
  // forward declaration
  class TimInt;

  /*!
  \brief abstract interface of time integration strategies for particle problems

  To keep the particle time integrator classes as plain as possible, several algorithmic parts focusing on the computation
  of kinematic and kinetic quantities have been encapsulated within separate time integration strategy classes. Thereby,
  a bunch of unhandy if-else selections within the time integrator classes themselves can be circumvented. Each specific
  time integration strategy class has to be implemented in a subclass derived from this abstract, purely virtual interface.

  \author fang (fang@lnm.mw.tum.de)

  \date 10/17
  */

  class TimIntStrategyBase
  {
    public:

      //! constructor
      explicit TimIntStrategyBase(
          TimInt* const   timint   //!< time integrator
          );

      //! destructor
      virtual ~TimIntStrategyBase(){return;};

      //! compute angular acceleration vector
      virtual void ComputeAngularAcceleration(
          Epetra_Vector&         angacc,     //!< global angular acceleration vector
          const Epetra_Vector&   m_contact   //!< global moment vector
          ) const = 0;

      //! compute inertia vector
      virtual void ComputeInertia() const = 0;

      //! compute kinetic energy
      virtual double ComputeKineticEnergy() const = 0;

      //! compute mass vector
      virtual void ComputeMass() const = 0;

      //! extract maximum particle radius
      virtual double MaxRadius() const = 0;

      //! output particle orientation
      virtual void OutputOrientation() const = 0;

      //! predict or correct angular velocity vector
      virtual void PredictOrCorrectAngularVelocity(
          const double           dt,          //!< time step size
          const Epetra_Vector&   m_contact,   //!< global moment vector
          const Epetra_Vector&   orient       //!< particle orientation vector at time t_n
          ) const = 0;

      //! update particle orientation vector
      virtual void RotateOrientVector(
          const double   dt   //!< time step size
          ) const = 0;

      //! set initial particle radii
      virtual void SetInitialRadii() const = 0;

      //! update map of inertia vector
      virtual void UpdateInertiaVectorMap() const = 0;

    protected:

      //! time integrator
      TimInt* const timint_;
  }; // class TimIntStrategyBase


  /*!
  \brief time integration strategy for spherical particles

  \author fang (fang@lnm.mw.tum.de)

  \date 10/17
  */

  class TimIntStrategySpheres : public TimIntStrategyBase
  {
    public:

      //! constructor
      explicit TimIntStrategySpheres(
          TimInt* const   timint   //!< time integrator
          );

      //! destructor
      virtual ~TimIntStrategySpheres(){return;};

      //! compute angular acceleration vector
      void ComputeAngularAcceleration(
          Epetra_Vector&         angacc,     //!< global angular acceleration vector
          const Epetra_Vector&   m_contact   //!< global moment vector
          ) const;

      //! compute inertia vector
      void ComputeInertia() const;

      //! compute kinetic energy
      double ComputeKineticEnergy() const;

      //! compute mass vector
      void ComputeMass() const;

      //! extract maximum particle radius
      double MaxRadius() const;

      //! output particle orientation
      void OutputOrientation() const;

      //! predict or correct angular velocity vector
      void PredictOrCorrectAngularVelocity(
          const double           dt,          //!< time step size
          const Epetra_Vector&   m_contact,   //!< global moment vector
          const Epetra_Vector&   orient       //!< particle orientation vector at time t_n
          ) const;

      //! update particle orientation vector
      void RotateOrientVector(
          const double   dt   //!< time step size
          ) const;

      //! set initial particle radii
      void SetInitialRadii() const;

      //! update map of inertia vector
      void UpdateInertiaVectorMap() const;
  }; // class TimIntStrategySpheres


  /*!
  \brief time integration strategy for ellipsoidal particles

  \author fang (fang@lnm.mw.tum.de)

  \date 10/17
  */

  class TimIntStrategyEllipsoids : public TimIntStrategyBase
  {
    public:

      //! constructor
      explicit TimIntStrategyEllipsoids(
          TimInt* const   timint   //!< time integrator
          );

      //! destructor
      virtual ~TimIntStrategyEllipsoids(){return;};

      //! compute angular acceleration vector
      void ComputeAngularAcceleration(
          Epetra_Vector&         angacc,     //!< global angular acceleration vector
          const Epetra_Vector&   m_contact   //!< global moment vector
          ) const;

      //! compute inertia vector
      void ComputeInertia() const;

      //! compute kinetic energy
      double ComputeKineticEnergy() const;

      //! compute mass vector
      void ComputeMass() const;

      //! extract maximum particle radius
      double MaxRadius() const;

      //! output particle orientation
      void OutputOrientation() const;

      //! predict or correct angular velocity vector
      void PredictOrCorrectAngularVelocity(
          const double           dt,          //!< time step size
          const Epetra_Vector&   m_contact,   //!< global moment vector
          const Epetra_Vector&   orient       //!< particle orientation vector at time t_n
          ) const;

      //! update particle orientation vector
      void RotateOrientVector(
          const double   dt   //!< time step size
          ) const;

      //! set initial particle radii
      void SetInitialRadii() const;

      //! update map of inertia vector
      void UpdateInertiaVectorMap() const;

    private:

      //! semi-axes of ellipsoidal particles
      LINALG::Matrix<3,1> semiaxes_;
  }; // class TimIntStrategyEllipsoids
} // namespace PARTICLE
#endif // #ifndef PARTICLE_TIMINT_STRATEGY_H
