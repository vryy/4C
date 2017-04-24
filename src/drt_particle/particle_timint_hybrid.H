/*----------------------------------------------------------------------*/
/*!
\file particle_timint_hybrid.H

\brief Hybrid particle time integration

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_HYBRID_H
#define PARTICLE_TIMINT_HYBRID_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Front-end for particle dynamics with hybrid time integrators
   *
   * <h3> About </h3>
   * This object bridges the gap between the base time integator PARITCLE::TimInt
   * and the specific implementation of hybrid time integrators.
   *
   * \author cattabiani
   * \date 01/17
 */
  class TimIntHybrid : public TimInt
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntHybrid
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntHybrid(){ ; }

    //@}

    /// update at the end of the time step
    virtual void Update()
    {
      UpdateStepState();
      UpdateStepTime();
      return;
    }

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;

    //@}

  };  // class TimIntHybrid

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_HYBRID_H_ */
