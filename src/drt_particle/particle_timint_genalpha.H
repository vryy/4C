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
      Teuchos::RCP<LINALG::Solver>& solver, //!< solver
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
    virtual void PredictMidState(const Teuchos::RCP<Epetra_Vector> disn);

    //! Evaluate new-state vectors by Newmark's formulae
    virtual void PredictNewState(const Teuchos::RCP<Epetra_Vector> disn, const bool disAreEqual = false);

    //! Evaluate the residual (acceleration)
    virtual void ResAcc();

    //! Evaluate the gradient of the residual (acceleration)
    virtual void GradResAcc();

    //! Compute deltaDis and update disn - Newton-Rhapson iteration
    virtual void CorrectDis(bool checkDis);

    //! Set up the additional state vectors required in the genAlpha algorithm
    virtual void SetupStateVectors();

    //! Compute and updates in the interactionHandler all the relevant dis-related state vectors
    virtual void ComputeAndSetDisRelatedStateVectors();

    //! debug too to check the gradResAcc against its residual difference definition
    virtual Teuchos::RCP<LINALG::SparseMatrix> FDGradResAcc();

    //! solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //@}

    //! mid-average type more at MidAverageEnum
    enum INPAR::PARTICLE::MidAverageEnum midavg_;

    //! Approximation type
    enum INPAR::PARTICLE::GAapproxType gradResAccApproxType_;

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
    const double tol_; //! gen-alpha residual tollerance
    const int maxIt_; //! gen-alpha maximum number of iterations




    Teuchos::RCP<Epetra_Vector> mGradW_;
    //! m * hess(weightFunction) for meshfree methods (not needed for now)
    //Teuchos::RCP<Epetra_MultiVector> mHessW_;
    //! Displacements \f${D}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> dism_;
    //! Velocities \f${V}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> velm_;
    //! Accelerations \f${A}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> accm_;
    //! Acceleration residue \f${R}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<Epetra_Vector> resAcc_;
    //! Gradient of the acceleration residue \f$\grad{R}_{n+1-\alpha f}\f$ - required for genalpha implicit scheme
    Teuchos::RCP<LINALG::SparseMatrix> gradResAcc_;

  };  // class TimIntGenAlpha

}  // namespace PARTICLE

#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_ */

