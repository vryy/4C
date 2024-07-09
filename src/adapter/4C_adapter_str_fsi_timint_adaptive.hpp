/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for time step size adaptivity within monolithic FSI

\level 2


*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_FSI_TIMINT_ADAPTIVE_HPP
#define FOUR_C_ADAPTER_STR_FSI_TIMINT_ADAPTIVE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_timint_adaptive.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations:
namespace Solid
{
  class TimInt;
  class TimAda;
}  // namespace Solid


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace Adapter
{
  /*====================================================================*/
  /*!
   * \brief Structure field adapter for time step size adaptivity within monolithic FSI
   *
   * Use this adapter in case you want to do monolithic FSI with time step size
   * adaptivity. By inheritance, we combine FSI functionalities with structural
   * time adaptivity. The FSI stuff is inherited from Adapter::FSIStructureWrapper
   * and the time adaptivity from Adapter::StructureTimIntAda
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
    StructureFSITimIntAda(Teuchos::RCP<Solid::TimAda> sta, Teuchos::RCP<Structure> sti);

    //! Do one time step with auxiliary time integration scheme
    virtual void time_step_auxiliar();

    //! Indicate norms of local discretization error
    virtual void indicate_error_norms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
    );

    //! Calculate time step size suggestion
    virtual double calculate_dt(const double norm);

    //! Get time step size of adaptive structural time integrator
    double dt() const override;

    //! Get target time \f$t_{n+1}\f$ of current time step
    double time() const override;

    //! Set new time step size
    void set_dt(const double dtnew) override;

    //! Update step size
    virtual void update_step_size(const double dtnew);

    //! Reset certain quantities to prepare repetition of current time step
    void reset_step() override;

    //! return pointer to structure time integration
    Teuchos::RCP<Structure> get_str_tim_int_ptr() { return str_time_integrator_; };

   private:
    //! Indicate local discretization error
    void indicate_errors(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
    );

    enum Inpar::Solid::VectorNorm errnorm_;  //!< norm for local error vector

    int numdbcdofs_;       ///< number of DOFs with Dirichlet boundary condition
    int numdbcfsidofs_;    ///< number of interface DOFs with Dirichlet boundary condition
    int numdbcinnerdofs_;  ///< number of inner DOFs with Dirichlet boundary condition

    Teuchos::RCP<Structure> str_time_integrator_;  ///< pointer to the structural time integrator

  };  // class StructureFSITimIntAda

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
