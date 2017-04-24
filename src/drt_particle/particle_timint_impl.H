/*----------------------------------------------------------------------*/
/*!
\file particle_timint_impl.H

\brief Implicit particle time integration

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_IMPL_H
#define PARTICLE_TIMINT_IMPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Front-end for particle dynamics with implicit time integrators
   *
   * <h3> About </h3>
   * This object bridges the gap between the base time integator PARITCLE::TimInt
   * and the specific implementation of implicit time integrators.
   *
   * \author cattabiani
   * \date 01/17
 */
  class TimIntImpl : public TimInt
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntImpl
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntImpl(){ ; }

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

  private :

    /// tolerance for the error in the iterative loop
    const double errToll_;

    /// max number of iterations of the iterative loop
    const double iterMax_;

  };  // class TimIntimpl

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_IMPL_H_ */
