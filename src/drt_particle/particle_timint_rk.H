/*----------------------------------------------------------------------*/
/*! \file
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

\level 3
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_RK_H
#define PARTICLE_TIMINT_RK_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{
  /*====================================================================*/
  /*!
   * \brief Runge-Kutta:         2nd/4th order accurate,
   *                             explicit time integrator
   * \author rasthofer
   * \date 01/14
   */
  class TimIntRK : public PARTICLE::TimIntExpl
  {
   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntRK(const Teuchos::ParameterList& ioparams,      //!< ioflags
        const Teuchos::ParameterList& particledynparams,  //!< input parameters
        const Teuchos::ParameterList& xparams,            //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,         //!< current discretisation
        Teuchos::RCP<IO::DiscretizationWriter> output     //!< the output
    );

    //! Destructor
    virtual ~TimIntRK() { ; }

    //@}

    //! @name Actions
    //@{

    //! do time integration of single time step
    virtual int IntegrateStep();

    //@}

   protected:
    //! Runge-Kutta scheme 2nd order
    void Integrate_RK_Second();

    //! Runge-Kutta scheme 4th order
    void Integrate_RK_Fourth();

   private:
    //! output displacement
    virtual void OutputDisplacement() const;

    //! @name General method parameters
    //@{

    INPAR::PARTICLEOLD::DynamicType rk_scheme_;

    //@}

  };  // class TimIntRK

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_RK_H
