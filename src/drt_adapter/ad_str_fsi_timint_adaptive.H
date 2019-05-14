/*----------------------------------------------------------------------*/
/*!

\brief Structure field adapter for time step size adaptivity within monolithic FSI

\level 2

\maintainer Matthias Mayr

*/

/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_FSI_TIMINT_ADAPT_H
#define ADAPTER_STRUCTURE_FSI_TIMINT_ADAPT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsiwrapper.H"
#include "ad_str_timint_adaptive.H"

// forward declarations:
namespace STR
{
  class TimInt;
  class TimAda;
}  // namespace STR


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER
{
  /*====================================================================*/
  /*!
   * \brief Structure field adapter for time step size adaptivity within monolithic FSI
   *
   * Use this adapter in case you want to do monolithic FSI with time step size
   * adaptivity. By inheritance, we combine FSI functionalities with structural
   * time adaptivity. The FSI stuff is inherited from ADAPTER::FSIStructureWrapper
   * and the time adaptivity from ADAPTER::StructureTimIntAda
   *
   * The time loop is implemented in FSI::Monolithic, which requires
   * error estimation and time step size calculation based on the structure field.
   * For error estimation and time step size calculation we want to use the standard
   * structural time adaptivity routines. Though, the decision, whether a time step
   * has to be repeated, has to be made by the FSI algorithm.
   *
   * \sa FSIStructureWrapper
   * \sa StructureTimIntAda
   *
   * \author mayr.mt
   * \date 12/2013
   */
  class StructureFSITimIntAda : virtual public FSIStructureWrapper,
                                virtual public StructureTimIntAda
  {
   public:
    //! Constructor
    StructureFSITimIntAda(Teuchos::RCP<STR::TimAda> sta, Teuchos::RCP<Structure> sti);

    //! Destructor
    virtual ~StructureFSITimIntAda() { ; }

    //! Do one time step with auxiliary time integration scheme
    virtual void TimeStepAuxiliar();

    //! Indicate norms of local discretization error
    virtual void IndicateErrorNorms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
    );

    //! Calculate time step size suggestion
    virtual double CalculateDt(const double norm);

    //! Get time step size of adaptive structural time integrator
    virtual double Dt() const;

    //! Get target time \f$t_{n+1}\f$ of current time step
    virtual double Time() const;

    //! Set new time step size
    virtual void SetDt(const double dtnew);

    //! Update step size
    virtual void UpdateStepSize(const double dtnew);

    //! Reset certain quantities to prepare repetition of current time step
    virtual void ResetStep();

    //! return pointer to structure time integration
    Teuchos::RCP<Structure> GetStrTimIntPtr() { return str_time_integrator_; };

   private:
    //! Indicate local discretization error
    void IndicateErrors(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
    );

    enum INPAR::STR::VectorNorm errnorm_;  //!< norm for local error vector

    int numdbcdofs_;       ///< number of DOFs with Dirichlet boundary condition
    int numdbcfsidofs_;    ///< number of interface DOFs with Dirichlet boundary condition
    int numdbcinnerdofs_;  ///< number of inner DOFs with Dirichlet boundary condition

    Teuchos::RCP<Structure> str_time_integrator_;  ///< pointer to the structural time integrator

  };  // class StructureFSITimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
