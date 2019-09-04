/*----------------------------------------------------------------------*/
/*! \file

\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

\level 2

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_CENTRDIFF_H
#define PARTICLE_TIMINT_CENTRDIFF_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{
  /*====================================================================*/
  /*!
   * \brief Central differences: 2nd order accurate,
   *                             explicit time integrator,
   *                             1-step method
   * \author hammerl
   * \date 09/12
   */
  class TimIntCentrDiff : public TimIntExpl
  {
   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntCentrDiff(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& particledynparams,     //!< input parameters
        const Teuchos::ParameterList& xparams,               //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,            //!< current discretisation
        Teuchos::RCP<IO::DiscretizationWriter> output        //!< the output
    );

    //! Destructor
    virtual ~TimIntCentrDiff() { ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Get extreme values of gap etc.
    virtual void GetExtremeValues();

    //@}

  };  // class TimIntCentrDiff

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_CENTRDIFF_H
