/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for partitioned poroelasticity algorithms

\level 2

 *-----------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_PARTITIONED_HPP
#define FOUR_C_POROELAST_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_inpar_poroelast.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_poroelast_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroElast
{
  //! base class of all monolithic Poroelasticity algorithms
  class Partitioned : public PoroBase
  {
   public:
    //! create using a Epetra_Comm
    explicit Partitioned(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter);

    //! proceed one time step (prepare, solve, update)
    void do_time_step() override;

    //! initialise system
    void SetupSystem() override;

    //! dof row map of Structure field
    Teuchos::RCP<const Epetra_Map> DofRowMapStructure() override;

    //! dof row map of Fluid field
    Teuchos::RCP<const Epetra_Map> DofRowMapFluid() override;

   protected:
    //! prepare new time step
    void prepare_time_step() override;

    //! solve one time step of structural problem
    void do_struct_step();

    //! solve one time step of fluid problem
    void do_fluid_step();

    //! solve one time step (iteration between fields)
    void Solve() override;

    //! update and write output to screen and files after solved time step
    void update_and_output();

    //! convergence check of outer loop
    bool convergence_check(int itnum);

    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;

    Teuchos::RCP<Epetra_Vector> fluidveln_;  //!< global fluid velocities and pressures
  };

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
