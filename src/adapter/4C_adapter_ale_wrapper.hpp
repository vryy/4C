/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 2

 */
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_ALE_WRAPPER_HPP
#define FOUR_C_ADAPTER_ALE_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_ale_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /*! \brief Just a wrapper that does nothing, meant to be derived from
   *
   *  This wrapper just encapsulated the Adapter::Ale and implements all
   *  routines that are pure virtual in Adapter::Ale. For a specific ALE adapter
   *  just derive from this one and overload those routines you need with your
   *  problem specific routine.
   *
   *  \author mayr.mt \date 10/2014
   */
  class AleWrapper : public Ale
  {
   public:
    //! constructor
    explicit AleWrapper(Teuchos::RCP<Ale> ale) : ale_(ale) {}

    //! @name Vector access
    //@{

    //! initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> initial_guess() const override
    {
      return ale_->initial_guess();
    }

    //! right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> rhs() const override { return ale_->rhs(); }

    //! unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> dispnp() const override { return ale_->dispnp(); }

    //! known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> dispn() const override { return ale_->dispn(); }

    //@}

    //! @name Misc
    //@{

    //! dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> dof_row_map() const override { return ale_->dof_row_map(); }

    //! direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return ale_->system_matrix();
    }

    //! direct access to system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
    {
      return ale_->block_system_matrix();
    }

    //! access to locsys manager
    Teuchos::RCP<Core::Conditions::LocsysManager> locsys_manager() override
    {
      return ale_->locsys_manager();
    }

    //! direct access to discretization
    Teuchos::RCP<const Core::FE::Discretization> discretization() const override
    {
      return ale_->discretization();
    }

    /// writing access to discretization
    Teuchos::RCP<Core::FE::Discretization> write_access_discretization() override
    {
      return ale_->write_access_discretization();
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor()
    {
      return ale_->get_dbc_map_extractor(ALE::UTILS::MapExtractor::dbc_set_std);
    }

    //! Return MapExtractor for Dirichlet boundary conditions in case of non-standard Dirichlet sets
    Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type  ///< type of dbc set
        ) override
    {
      return ale_->get_dbc_map_extractor(dbc_type);
    }

    //! reset state vectors to zero
    void reset() override { ale_->reset(); }

    //! reset last time step, needed for time step size adaptivity of FSI
    void reset_step() override { ale_->reset_step(); }

    //@}

    //! @name Time step helpers
    //@{

    void reset_time(const double dtold) override { ale_->reset_time(dtold); }
    //! Return target time \f$t_{n+1}\f$
    double time() const override { return ale_->time(); }

    //! Return target step counter \f$step_{n+1}\f$
    double step() const override { return ale_->step(); }

    //! get time step size \f$\Delta t_n\f$
    double dt() const override { return ale_->dt(); }

    //! integrate from t1 to t2
    int integrate() override { return ale_->integrate(); }

    void time_step(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
                       ALE::UTILS::MapExtractor::dbc_set_std) override
    {
      ale_->time_step(dbc_type);
      return;
    }

    //! set time step size
    void set_dt(const double dtnew  ///< new time step size (to be set)
        ) override
    {
      ale_->set_dt(dtnew);
    }

    //! Set time and step
    void set_time_step(const double time,  ///< simulation time (to be set)
        const int step                     ///< step number (to be set)
        ) override
    {
      ale_->set_time_step(time, step);
    }

    //! start new time step
    void prepare_time_step() override { ale_->prepare_time_step(); }

    //! update displacement and evaluate elements
    virtual void evaluate(Teuchos::RCP<const Epetra_Vector> stepinc =
                              Teuchos::null  ///< step increment such that \f$ x_{n+1}^{k+1} =
                                             ///< x_{n}^{converged}+ stepinc \f$
    )
    {
      evaluate(stepinc, ALE::UTILS::MapExtractor::dbc_set_std);
    }

    //! update displacement and evaluate elements
    void evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc,  ///< step increment such that \f$ x_{n+1}^{k+1}
                                                    ///< = x_{n}^{converged}+ stepinc \f$
        ALE::UTILS::MapExtractor::AleDBCSetType
            dbc_type  ///< application-specific type of Dirichlet set
        ) override
    {
      ale_->evaluate(stepinc, dbc_type);
    }

    //! update at time step end
    void update() override { ale_->update(); }

    //! update at time step end
    void update_iter() override { ale_->update_iter(); }

    //! output results
    void output() override { return ale_->output(); }


    //! read restart information for given time step \p step
    void read_restart(const int step  ///< step number to read restart from
        ) override
    {
      return ale_->read_restart(step);
    }

    //@}

    /// setup Dirichlet boundary condition map extractor
    void setup_dbc_map_ex(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std,  //!< application-specific type of Dirichlet set
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null,  //!< interface for creation of additional, application-specific
                            //!< Dirichlet map extractors
        Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface =
            Teuchos::null  //!< interface for creation of a Dirichlet map extractor, taylored to
                           //!< XFFSI
        ) override
    {
      ale_->setup_dbc_map_ex(dbc_type, interface, xff_interface);
    }

    //! @name Solver calls
    //@{

    //! nonlinear solve
    int solve() override { return ale_->solve(); }

    //! Access to linear solver of ALE field
    Teuchos::RCP<Core::LinAlg::Solver> linear_solver() override { return ale_->linear_solver(); }

    //@}

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    //! write access to extract displacements at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> write_access_dispnp() const override
    {
      return ale_->write_access_dispnp();
    }

    //@}

    //! create result test for encapsulated structure algorithm
    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override
    {
      return ale_->create_field_test();
    }

    /*! \brief Create Systemmatrix
     *
     * We allocate the Core::LINALG object just once, the result is an empty Core::LINALG object.
     * Evaluate has to be called separately.
     *
     */
    void create_system_matrix(
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null  ///< Blocksparsematrix if an interface is passed, SparseMatrix otherwise
        ) override
    {
      ale_->create_system_matrix(interface);
    }

    //! update slave dofs for fsi simulations with ale mesh tying
    void update_slave_dof(Teuchos::RCP<Epetra_Vector>& a) override { ale_->update_slave_dof(a); }

   private:
    Teuchos::RCP<Ale> ale_;  //!< underlying ALE time integration
  };

  //! Calculate increments from absolute values
  class AleNOXCorrectionWrapper : public AleWrapper  // ToDo (mayr) Do we really need this?
  {
   public:
    explicit AleNOXCorrectionWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale) {}

    //! Prepare time step
    void prepare_time_step() override;

    /*! \brief evaluate() routine that can handle NOX step increments
     *
     *  We deal with NOX step increments by computing the last iteration increment
     *  needed for the ALE evaluate() call.
     *
     *  The field solver always expects an iteration increment only. And
     *  there are Dirichlet conditions that need to be preserved. So take
     *  the sum of increments we get from NOX and apply the latest iteration
     *  increment only.
     *  Naming:
     *
     *  \f$x^{n+1}_{i+1} = x^{n+1}_i + iterinc\f$  (sometimes referred to as residual increment),
     * and
     *
     *  \f$x^{n+1}_{i+1} = x^n + stepinc\f$
     *
     *  \author mayr.mt \date 10/2014
     */
    void evaluate(Teuchos::RCP<const Epetra_Vector> stepinc  ///< step increment
        ) override;

   private:
    //! sum of displacement increments already applied,
    //!
    //! there are two increments around
    //!
    //! x^n+1_i+1 = x^n+1_i + stepinc  (also referred to as residual increment)
    //!
    //! x^n+1_i+1 = x^n     + disstepinc
    Teuchos::RCP<Epetra_Vector> stepinc_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
