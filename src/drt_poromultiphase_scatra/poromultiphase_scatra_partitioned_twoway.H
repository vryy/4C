/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_partitioned_twoway.H

 \brief two-way coupled partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_


#include "poromultiphase_scatra_partitioned.H"

namespace POROMULTIPHASESCATRA
{
  //! Base class of all partitioned solid-scatra algorithms --> virtual
  class PoroMultiPhaseScaTraPartitionedTwoWay : public PoroMultiPhaseScaTraPartitioned
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitionedTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraPartitionedTwoWay(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra);

    /// setup
    virtual void SetupSystem();

    /// setup solver (only needed for poromultiphase monolithic coupling)
    virtual void SetupSolver();

    /// time step of coupled problem
    void TimeStep() { return Solve(); };

    /// print header
    void PrintHeaderPartitioned();

    /// print header
    void IterUpdateStates();

    //! perform iteration loop between fields
    virtual void Solve() = 0;


   protected:
    //! perform iteration step of structure-fluid field
    void DoPoroStep();

    //! perform iteration step of scatra field
    void DoScatraStep();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;
    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;


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

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraPartitionedTwoWayNested(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra);

    //! perform iteration loop between fields
    virtual void Solve();

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

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraPartitionedTwoWaySequential(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra);

    //! perform iteration loop between fields
    virtual void Solve();

  };  // PoroMultiPhaseScatraPartitionedTwoWayNested


}  // namespace POROMULTIPHASESCATRA



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_ */
