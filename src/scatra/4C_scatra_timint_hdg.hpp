/*----------------------------------------------------------------------*/
/*! \file
\brief Time integration class for HDG discretization for scalar transport

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_HDG_HPP
#define FOUR_C_SCATRA_TIMINT_HDG_HPP

#include "4C_config.hpp"

#include "4C_lib_discret_hdg.hpp"
#include "4C_scatra_timint_genalpha.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  /*!
  \brief time integration for HDG scatra
  */
  class TimIntHDG : public TimIntGenAlpha
  {
   public:
    //! standard constructor
    TimIntHDG(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Teuchos::ParameterList>& extraparams,
        Teuchos::RCP<CORE::IO::DiscretizationWriter> output);

    //! setup
    void Setup() override;

    //! set theta_ to its value, dependent on integration method for GenAlpha and BDF2
    virtual void SetTheta();

    //! set states in the time integration schemes: additional vectors for HDG
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    //! set the part of the right hand side belonging to the last time step
    void set_old_part_of_righthandside() override;

    //! update the solution after convergence of the nonlinear iteration,
    //! current solution becomes old solution of next time step
    void Update() override;

    //  //! Initialization procedure before the first time step is done
    //  void prepare_first_time_step ();

    //! update configuration and output to file/screen
    void output_state() override;

    void write_restart() const override;

    //! read restart
    void read_restart(
        const int step, Teuchos::RCP<CORE::IO::InputControl> input = Teuchos::null) override;

    //! set the initial scalar field phi
    void SetInitialField(const INPAR::SCATRA::InitialField init,  //!< type of initial field
        const int startfuncno                                     //!< number of spatial function
        ) override;

    //! accessor to interior concentrations
    virtual Teuchos::RCP<Epetra_Vector> ReturnIntPhinp() { return intphinp_; }
    virtual Teuchos::RCP<Epetra_Vector> ReturnIntPhin() { return intphin_; }
    //  virtual Teuchos::RCP<Epetra_Vector>  ReturnIntPhinm(){return intphinm_;}

    virtual Teuchos::RCP<Epetra_Vector> InterpolatedPhinp() const { return interpolatedPhinp_; }

    //! prepare time loop
    void prepare_time_loop() override;

    /*!
    \brief Compare the numerical solution to the analytical one.
    */
    virtual Teuchos::RCP<CORE::LINALG::SerialDenseVector> compute_error() const;

    Teuchos::RCP<CORE::UTILS::ResultTest> create_sca_tra_field_test() override;

   protected:
    //! copy constructor
    TimIntHDG(const TimIntHDG& old);

    //! update time derivative for generalized-alpha time integration
    virtual void gen_alpha_compute_time_derivative();

    //! compute values at intermediate time steps for gen.-alpha
    virtual void gen_alpha_intermediate_values();

    //! number of dofset for interior variables
    int nds_intvar_;

    //! @name concentration and concentration gradient at different times for element interior for
    //! HDG
    //@{
    Teuchos::RCP<Epetra_Vector> intphinp_;  //!< concentration at time \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> intphin_;   //!< concentration at time \f$t^{n}\f$
    //  Teuchos::RCP<Epetra_Vector> intphinm_;   //!< concentration at time \f$t^{n-1}\f$
    //  Teuchos::RCP<Epetra_Vector> intphiaf_;   //!< concentration at time \f$t^{n+\alpha_F}\f$
    //  Teuchos::RCP<Epetra_Vector> intphiam_;   //!< concentration at time \f$t^{n+\alpha_M}\f$
    //@}

    //! @name scalar time derivative of concentration and concentration gradient
    //! at time n+1, n and n+alpha_M/(n+alpha_M/n) and n-1 for element interior in HDG
    //@{
    //  Teuchos::RCP<Epetra_Vector> intphidtnp_;   //!< time derivative at time \f$t^{n+1}\f$
    //  Teuchos::RCP<Epetra_Vector> intphidtn_;    //!< time derivative at time \f$t^{n}\f$
    //  Teuchos::RCP<Epetra_Vector> intphidtnm_;   //!< time derivative at time \f$t^{n-1}\f$
    //  Teuchos::RCP<Epetra_Vector> intphidtam_;   //!< time derivative at time \f$t^{n+\alpha_M}\f$
    //@}

    //! @name other HDG-specific auxiliary vectors
    //@{
    Teuchos::RCP<Epetra_Vector>
        interpolatedPhinp_;  //!< concentrations for output at time \f$t^{n+1}\f$
    //@}

    //! calculate intermediate solution
    void compute_intermediate_values() override;

    //! compute values at the interior of the elements
    void compute_interior_values() override;

    //! update interior variables
    virtual void update_interior_variables(Teuchos::RCP<Epetra_Vector> updatevector);

    //! write problem specific output
    virtual void write_problem_specific_output(Teuchos::RCP<Epetra_Vector> interpolatedPhi)
    {
      return;
    }

    //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
    void calc_initial_time_derivative() override { return; };

    //! fd check
    void fd_check() override;

    //! calculation of error with reference to analytical solution during the simulation
    void evaluate_error_compared_to_analytical_sol() override;

    //! adapt degree of test functions and change dofsets accordingly
    virtual void adapt_degree();

    //! adapt variable vectors required due to the change of the degrees of the test functions
    virtual void adapt_variable_vector(Teuchos::RCP<Epetra_Vector> phi_new,
        Teuchos::RCP<Epetra_Vector> phi_old, Teuchos::RCP<Epetra_Vector> intphi_new,
        Teuchos::RCP<Epetra_Vector> intphi_old, int nds_var_old, int nds_intvar_old,
        std::vector<CORE::Elements::Element::LocationArray> la_old);

    //! calculate matrices on element
    virtual void calc_mat_initial();

    //! chooses the assemble process (assemble matrix and rhs or only rhs)
    void assemble_mat_and_rhs() override;

    //! contains the assembly process only for rhs
    void assemble_rhs();

    //! pack material
    virtual void pack_material() { return; };

    //! adapt material
    virtual void unpack_material() { return; };

    //! project material field
    virtual void project_material() { return; };

   private:
    //! time algorithm flag actually set (we internally reset it)
    INPAR::SCATRA::TimeIntegrationScheme timealgoset_;

    //! @name time stepping variable
    bool startalgo_;  //!< flag for starting algorithm

    //! time-integration-scheme factor theta
    double theta_;

    //! activation_time at times n+1
    Teuchos::RCP<Epetra_Vector> activation_time_interpol_np_;

    //! HDG discretization
    DRT::DiscretizationHDG* hdgdis_;

    //! p-adativitity
    bool padaptivity_;

    //! error tolerance for p-adativitity
    double padapterrortol_;

    //! error tolerance base for p-adativitity (delta p=log_padapterrorbase_(E/padapterrortol_))
    double padapterrorbase_;


    //! max. degree of shape functions for p-adativitity
    double padaptdegreemax_;

    //! element degree
    Teuchos::RCP<Epetra_Vector> elementdegree_;

  };  // class TimIntHDG
}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
