/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_HPP


#include "4C_config.hpp"

#include "4C_poromultiphase_scatra_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroMultiPhaseScaTra
{
  //! Base class of all partitioned solid-scatra algorithms --> virtual
  class PoroMultiPhaseScaTraPartitionedTwoWay : public PoroMultiPhaseScaTraPartitioned
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitionedTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    /// setup
    void setup_system() override;

    /// setup solver (only needed for poromultiphase monolithic coupling)
    void setup_solver() override;

    /// time step of coupled problem
    void time_step() override { return solve(); };

    /// print header
    void print_header_partitioned();

    /// print header
    void iter_update_states();

    //! perform iteration loop between fields
    virtual void solve() = 0;


   protected:
    //! perform iteration step of structure-fluid field
    void do_poro_step();

    //! perform iteration step of scatra field
    void do_scatra_step();

    //! convergence check of outer loop
    bool convergence_check(int itnum);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;
    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;
    //! artery scatra increment of the outer loop
    Teuchos::RCP<Epetra_Vector> artscaincnp_;
    //! artery pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> arterypressincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
    //! is artery coupling active
    bool artery_coupling_active_;


  };  // PoroMultiPhaseScatraPartitionedTwoWay

  //! Nested partitioned solution algorithm
  // +--------------------------+           +----------+
  // |         ---->            |  ------>  |          |
  // |  fluid        structure  |           | ScaTra   |
  // |         <----            |  <------  |          |
  // +--------------------------+           +----------+

  class PoroMultiPhaseScaTraPartitionedTwoWayNested : public PoroMultiPhaseScaTraPartitionedTwoWay
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitionedTwoWayNested(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    //! perform iteration loop between fields
    void solve() override;

  };  // PoroMultiPhaseScatraPartitionedTwoWayNested

  //! Sequential partitioned solution algorithm
  // +-----------+          +-----------+           +-----------+
  // |           |  ----->  |           | --------> |           |
  // |   fluid   |          | structure |           |  ScaTra   |
  // |           |          |           |           |           |
  // +-----------+          +-----------+           +-----------+
  //      ^                                              |
  //      |----------------------------------------------+

  class PoroMultiPhaseScaTraPartitionedTwoWaySequential
      : public PoroMultiPhaseScaTraPartitionedTwoWay
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitionedTwoWaySequential(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    //! perform iteration loop between fields
    void solve() override;

  };  // PoroMultiPhaseScatraPartitionedTwoWayNested


}  // namespace PoroMultiPhaseScaTra



FOUR_C_NAMESPACE_CLOSE

#endif
