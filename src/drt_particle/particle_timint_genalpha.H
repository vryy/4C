/*----------------------------------------------------------------------*/
/*!
\file particle_timint_genalpha.H

\brief Implicit particle time integration scheme

\level 2

<pre>
\maintainer Alessandro Cattabiani

</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_
#define SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_impl.H"
#include "../drt_inpar/inpar_particle.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Implicit time integration scheme (backward Euler?)
   *
   * \author Cattabiani
   */

  class TimIntGenAlpha : public TimIntImpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntGenAlpha
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntGenAlpha(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! overload to determine the initial accelerations
    virtual void DetermineMassDampConsistAccel();

    //! State vectors are updated according to the new distribution of particles (overload)
    virtual void UpdateStatesAfterParticleTransfer();

    //! Read and set restart state
    virtual void ReadRestartState();

    //@}

  protected:

    //! verify if given coefficients are in admissible range;
    //! prints also info to STDOUT
    virtual void VerifyCoeff();

    //! calculate coefficients from given spectral radius
    virtual void CalcCoeff();

    //! Evaluate mid-state vectors by averaging end-point vectors
    virtual void MidState();

    //! Evaluate new-state vectors by Newmark's formulae
    virtual void NewState(const bool disAreEqual = false);

    //! Evaluate the residual (acceleration)
    virtual void ResidualAcc();

    //! Set up the additional state vectors required in the genAlpha algorithm
    virtual void SetupStateVectors();

    //! mid-average type more at MidAverageEnum
    enum INPAR::PARTICLE::MidAverageEnum midavg_;
    //@}

    //! @name Key coefficients
    //! Please note, to obtain a second-order accurate scheme, you need
    //! to follow the following formulas in which \f$\rho_\infty\f$ is the
    //! spectral radius.
    //! \f[ \alpha_m = (2*\rho_\infty - 1)/(\rho_\infty + 1) \f]
    //! \f[ \alpha_f = \rho_\infty/(\rho_\infty + 1) \f]
    //! \f[ \beta = 1/4*(1 - \alpha_m + \alpha_f)^2 \mbox{(max. damp. of high-freq. modes)} \f]
    //! \f[ \gamma = 1/2 - \alpha_m + \alpha_f \f]
    //! The spectral radius is responsible for the magnitude of
    //! numerical dissipation introduced.
    //! For instance
    //! Without numerical dissipation at \f$\rho_\infty=1\f$
    //! \f[ \beta=0.25, \gamma=0.5, \alpha_m=0.5, \alpha_f=0.5 \f]
    //! Medium dissipation at \f$\rho_\infty=0.8\f$
    //! \f[ \beta=25/81, \gamma=11/18, \alpha_m=1/3, \alpha_f=4/9  \f]
    //! Strong numerical dissipation at \f$\rho_\infty=0.5\f$
    //! \f[ \beta=4/9, \gamma=10/12, \alpha_m=0, \alpha_f=1/3 \f]
    //@{
    double beta_;  //!< factor (0,1/2]
    double gamma_;  //!< factor (0,1]
    double alphaf_;  //!< factor [0,1)
    double alpham_;  //!< factor [-1,1)
    double rho_inf_; //!< factor[0,1]

    Teuchos::RCP<Epetra_Vector> mGradW_;
    //! m * hess(weightFunction) for meshfree methods
    Teuchos::RCP<Epetra_MultiVector> mHessW_;
    //! Displacements \f${D}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> dism_;
    //! Velocities \f${V}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> velm_;
    //! Accelerations \f${A}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> accm_;
    //! Acceleration residue \f${R}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> resAcc_;

  };  // class TimIntGenAlpha

}  // namespace PARTICLE

#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_ */

