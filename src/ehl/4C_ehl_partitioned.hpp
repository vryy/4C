/*--------------------------------------------------------------------------*/
/*! \file

\brief class for partitioned elastohydrodynamic lubrication (lubrication structure interaction)

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_EHL_PARTITIONED_HPP
#define FOUR_C_EHL_PARTITIONED_HPP


#include "4C_config.hpp"

#include "4C_ehl_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EHL
{
  class Partitioned : public Base
  {
   public:
    /// setup EHL algorithm
    explicit Partitioned(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string lubrication_disname);  // Problem builder

    /// setup of single fields (if needed)
    void SetupSystem() override{};

    /// time loop of coupled problem
    void Timeloop() override;

   protected:
    /// prepare time step of single fields
    void prepare_time_step() override;

    //! perform iteration loop between fields
    void outer_loop();

    //! update time step and print to screen
    void UpdateAndOutput();

    //! convergence check of outer loop
    bool convergence_check(int itnum);

    /// do one inner iteration loop step of the structure
    void do_struct_step();

    /// do one inner iteration loop step of the lubrication
    void DoLubricationStep();

    //! pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> preincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
  };
}  // namespace EHL

FOUR_C_NAMESPACE_CLOSE

#endif
