/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.H

\brief Particle time integration with central difference scheme 2nd order (explicit),
              also known as Velocity-Verlet algorithm

\level 2

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
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
    TimIntCentrDiff
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntCentrDiff(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! small support function for ComputeThermodynamics - update density
    void DensityUpdater(const double &thermalExpansion, const double &delta, double &density) {density /= (EffExpCoeff(thermalExpansion, delta)); }

    //! small support function for ComputeThermodynamics - update radius
    void RadiusUpdater(const double &thermalExpansion, const double &delta, double &radius) {radius *= std::pow(EffExpCoeff(thermalExpansion, delta), 1/3.0); }

    //@}

  protected:

    //! integrate step - displacement-related physics
    virtual void ComputeDisplacements();

    //! integrate step - thermodynamic-related physics
    virtual void ComputeThermodynamics1();


    //! integrate step - temperature related physics
    virtual int ComputeThermodynamics();

    //! expansion coefficient for radius and density (in case of thermodinamics)
    double EffExpCoeff(const double &thermalExpansion, const double &delta) {return thermalExpansion * delta + 1; }

  };  // class TimIntCentrDiff

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_CENTRDIFF_H
