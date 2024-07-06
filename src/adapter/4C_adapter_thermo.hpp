/*----------------------------------------------------------------------*/
/*! \file

\brief Thermo field adapter
\level 1



*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_THERMO_HPP
#define FOUR_C_ADAPTER_THERMO_HPP


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_thermo.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                      dano 02/11 |
 *----------------------------------------------------------------------*/
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
}  // namespace Core::LinAlg

namespace Core::Conditions
{
  class LocsysManager;
}


namespace CONTACT
{
  class NitscheStrategyTsi;
  class ParamsInterface;
}  // namespace CONTACT

namespace Adapter
{
  /// general thermal field interface
  /*!
  The point is to keep T(F)SI as far apart from our field solvers as
  possible. Each thermal field solver we want to use should get its own
  subclass of this. The T(F)SI algorithm should be able to extract all the
  information from the thermal field it needs using this interface.

  All T(F)SI algorithms use this adapter to communicate with the thermal
  field. There are different ways to use this adapter.

  In all cases you need to tell the thermal algorithm about your time
  step. Therefore prepare_time_step(), update() and output() must be called at
  the appropriate position in the TSI algorithm.

  <h3>Dirichlet-Neumann coupled TSI</h3>

  Dirichlet-Neumann coupled TSI will need to Solve() the linear thermal problem
  for each time step after the structure displacements/velocities have been
  applied (ApplyStructVariables()). Solve() will be called many times for each
  time step until the equilibrium is reached. The thermal algorithm has to
  preserve its state until update() is called.

  After each Solve() you get the new temperatures by Tempnp().

  <h3>Monolithic TSI</h3>

  Monolithic TSI is based on evaluate() of elements. This results in a new
  RHS() and a new SysMat(). Together with the initial_guess() these form the
  building blocks for a block based Newton's method.

  \warning Further cleanup is still needed.

  \sa Fluid, Ale, Structure
  \author cd, bborn
  \date 07/09
  */
  class Thermo
  {
   public:
    /// virtual to get polymorph destruction
    virtual ~Thermo() = default;

    /// @name Vector access
    //@{

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> initial_guess() = 0;

    /// RHS of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> rhs() = 0;

    /// unknown temperatures at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> tempnp() = 0;

    /// unknown temperatures at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> tempn() = 0;

    //@}

    /// @name Misc
    //@{

    /// DOF map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> dof_row_map() = 0;

    /// DOF map of vector of unknowns for multiple dofsets
    virtual Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) = 0;

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& domain_map() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<Core::FE::Discretization> discretization() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() = 0;

    //@}

    //! @name Time step helpers
    //@{
    //! Return current time \f$t_{n}\f$
    virtual double time_old() const = 0;

    //! Return target time \f$t_{n+1}\f$
    virtual double time() const = 0;

    /// Get upper limit of time range of interest
    virtual double get_time_end() const = 0;

    /// Get time step size \f$\Delta t_n\f$
    virtual double dt() const = 0;

    /// Return current step number $n$
    virtual int step_old() const = 0;

    /// Return current step number $n+1$
    virtual int step() const = 0;

    /// Get number of time steps
    virtual int num_step() const = 0;

    /// Set time step size for the current step
    virtual void set_dt(double timestepsize) = 0;

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    virtual void set_timen(const double time) = 0;

    /// Take the time and integrate (time loop)
    void integrate();

    /// tests if there are more time steps to do
    virtual bool not_finished() const = 0;

    /// start new time step
    virtual void prepare_time_step() = 0;

    /// evaluate residual at given temperature increment
    virtual void evaluate(Teuchos::RCP<const Epetra_Vector> tempi) = 0;

    /// evaluate residual at given temperature increment
    virtual void evaluate() = 0;

    /// update temperature increment after Newton step
    virtual void update_newton(Teuchos::RCP<const Epetra_Vector> tempi) = 0;

    /// update at time step end
    virtual void update() = 0;

    /// print info about finished time step
    virtual void print_step() = 0;

    //! Access to output object
    virtual Teuchos::RCP<Core::IO::DiscretizationWriter> disc_writer() = 0;

    /// prepare output
    virtual void prepare_output() = 0;

    /// output results
    virtual void output(bool forced_writerestart = false) = 0;

    /// read restart information for given time step
    virtual void read_restart(const int step) = 0;

    /// reset everything to beginning of time step, for adaptivity
    virtual void reset_step() = 0;

    //! store an RCP to the contact strategy constructed in the structural time integration
    virtual void set_nitsche_contact_strategy(
        Teuchos::RCP<CONTACT::NitscheStrategyTsi> strategy) = 0;

    //! store an RCP to the contact interface parameters in the structural time integration
    virtual void set_nitsche_contact_parameters(Teuchos::RCP<CONTACT::ParamsInterface> params) = 0;

    /// apply interface loads on the thermal field
    virtual void set_force_interface(Teuchos::RCP<Epetra_Vector> ithermoload) = 0;


    //@}

    //! @name Solver calls
    //@{

    /**
     * @brief Non-linear solve
     *
     * Do the nonlinear solve, i.e. (multiple) corrector for the time step.
     * All boundary conditions have been set.
     *
     * @return status of the solve, which can be used for adaptivity
     */
    virtual Inpar::THR::ConvergenceStatus solve() = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<Core::LinAlg::Solver> linear_solver() = 0;

    //@}

    //! @name Extract temperature values needed for TSI
    //@{
    /// extract temperatures for inserting in structure field
    virtual Teuchos::RCP<Epetra_Vector> write_access_tempn() = 0;

    /// extract current temperatures for inserting in structure field
    virtual Teuchos::RCP<Epetra_Vector> write_access_tempnp() = 0;
    //@}

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    virtual void prepare_partition_step() = 0;

    /// create result test for encapulated thermo algorithm
    virtual Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() = 0;
  };


  /// thermo field solver
  class ThermoBaseAlgorithm
  {
   public:
    /// constructor
    explicit ThermoBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis);

    /// virtual destructor to support polymorph destruction
    virtual ~ThermoBaseAlgorithm() = default;

    /// thermal field solver
    Thermo& thermo_field() { return *thermo_; }
    /// const version of thermal field solver
    const Thermo& thermo_field() const { return *thermo_; }
    /// Teuchos::rcp version of thermal field solver
    Teuchos::RCP<Thermo> thermo_fieldrcp() { return thermo_; }

   private:
    /// setup thermo algorithm
    void setup_thermo(
        const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis);

    /// setup thermo algorithm of THR::TimIntImpl type
    void setup_tim_int(const Teuchos::ParameterList& prbdyn, Inpar::THR::DynamicType timinttype,
        Teuchos::RCP<Core::FE::Discretization> actdis);

    /// thermal field solver
    Teuchos::RCP<Thermo> thermo_;

  };  // class ThermoBaseAlgorithm

}  // namespace Adapter


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
