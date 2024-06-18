/*----------------------------------------------------------------------*/
/*! \file

\brief  Connecting time-integration schemes for HDG with
        cardiac-monodomain-specific implementation

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_HPP
#define FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_cardiac_monodomain.hpp"
#include "4C_scatra_timint_hdg.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class TimIntCardiacMonodomainHDG : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntHDG
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainHDG(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! write current state to BINIO
    void output_state() override;

   protected:
    void element_material_time_update() override;

    //! write problem specific output
    void write_problem_specific_output(Teuchos::RCP<Epetra_Vector> interpolatedPhi) override;

    //! adapt material
    void pack_material() override;

    //! adapt material
    void unpack_material() override;

    //! project material field
    void project_material() override;

    //! read restart
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

   private:
    //! activation time
    Teuchos::RCP<Epetra_Vector> activation_time_interpol_;

    //! element data
    Teuchos::RCP<std::vector<char>> data_;


  };  // class TimIntCardiacMonodomainHDG

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
