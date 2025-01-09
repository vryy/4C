// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LUBRICATION_TIMINT_IMPLICIT_HPP
#define FOUR_C_LUBRICATION_TIMINT_IMPLICIT_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_lubrication.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Map.h>

FOUR_C_NAMESPACE_OPEN

/*==========================================================================*/
// Style guide                                                    nis Mar12
/*==========================================================================*/

/*--- set, prepare, and predict ------------------------------------------*/

/*--- calculate and update -----------------------------------------------*/

/*--- query and output ---------------------------------------------------*/



/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace Core::DOFSets
{
  class DofSet;
}  // namespace Core::DOFSets

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
  class KrylovProjector;
}  // namespace Core::LinAlg

namespace FLD
{
  class DynSmagFilter;
  class Vreman;
}  // namespace FLD

namespace Lubrication
{
  /*!
   * \brief implicit time integration for lubrication problems
   */

  class TimIntImpl
  {
   public:
    virtual std::shared_ptr<Core::IO::DiscretizationWriter> disc_writer() { return output_; }

    std::shared_ptr<Core::LinAlg::Vector<double>>& inf_gap_toggle() { return inf_gap_toggle_lub_; }

    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntImpl(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntImpl() = default;

    //! initialize time integration
    virtual void init();

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! set the nodal film height
    void set_height_field_pure_lub(const int nds);
    //! set the nodal film height
    void set_height_field(const int nds, std::shared_ptr<const Core::LinAlg::Vector<double>> gap);

    //! set the time derivative of the height (film thickness) by OST
    void set_height_dot_field(
        const int nds, std::shared_ptr<const Core::LinAlg::Vector<double>> heightdot);

    //! set relative tangential interface velocity for Reynolds equation
    void set_average_velocity_field_pure_lub(const int nds);
    void set_relative_velocity_field(
        const int nds, std::shared_ptr<const Core::LinAlg::Vector<double>> rel_vel);

    //! set average tangential interface velocity for Reynolds equation
    void set_average_velocity_field(
        const int nds, std::shared_ptr<const Core::LinAlg::Vector<double>> av_vel);

    //! add global state vectors specific for time-integration scheme
    virtual void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) = 0;

    //! prepare time loop
    virtual void prepare_time_loop();

    //! setup the variables to do a new time step
    virtual void prepare_time_step();

    //! initialization procedure prior to evaluation of first time step
    virtual void prepare_first_time_step();

    //! read restart data
    virtual void read_restart(int step) = 0;

    /*--- calculate and update -----------------------------------------------*/

    //! do time integration (time loop)
    virtual void time_loop();

    //! general solver call for coupled algorithms (decides if linear/nonlinear internally)
    virtual void solve();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void update(const int num = 0  //!< field number
        ) = 0;

    //! apply moving mesh data
    void apply_mesh_movement(
        std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp,  //!< displacement vector
        int nds  //!< number of the dofset the displacement state belongs to
    );

    //! calculate error compared to analytical solution
    virtual void evaluate_error_compared_to_analytical_sol();

    /*--- query and output ---------------------------------------------------*/

    //! print information about current time step to screen
    virtual void print_time_step_info();

    //! return system matrix downcasted as sparse matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix();

    //! update Newton step
    virtual void update_newton(std::shared_ptr<const Core::LinAlg::Vector<double>> prei);

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #prenp_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #predict_const_temp_consist_rate)
    virtual void update_iter_incrementally() = 0;

    //! Update iteration incrementally with prescribed residual
    //! pressures
    void update_iter_incrementally(const std::shared_ptr<const Core::LinAlg::Vector<double>>
            prei  //!< input residual pressures
    );

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic EHL accesses the linearised lubrication problem
    void evaluate();

    //! non-overlapping DOF map for multiple dofsets
    std::shared_ptr<const Epetra_Map> dof_row_map(unsigned nds = 0)
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map(nds);
      return std::make_shared<Epetra_Map>(*dofrowmap);
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() const
    {
      return dbcmaps_;
    }

    //! right-hand side alias the dynamic force residual
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() { return residual_; }

    //! return flag indicating if an incremental solution approach is used
    bool is_incremental() { return incremental_; }

    //! return discretization
    std::shared_ptr<Core::FE::Discretization> discretization() { return discret_; }

    //! output solution and restart data to file
    virtual void output(const int num = 0);

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- query and output ---------------------------------------------------*/

    //! return current time value
    double time() const { return time_; }

    //! return current step number
    int step() const { return step_; }

    //! return number of newton iterations in last timestep
    double iter_num() const { return iternum_; }

    //! return time step size
    double dt() const { return dta_; }

    /*========================================================================*/
    //! @name pressure degrees of freedom and related
    /*========================================================================*/

    /*--- query and output ---------------------------------------------------*/

    //! return pressure field pre at time n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> prenp() { return prenp_; }

    //! output mean values of pressure(s)
    virtual void output_mean_pressures(const int num = 0);

    //! output domain or boundary integrals, i.e., surface areas or volumes of specified nodesets
    void output_domain_or_boundary_integrals(const std::string condstring);

   protected:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! don't want = operator
    // TimIntImpl operator = (const TimIntImpl& old);

    //! don't want copy constructor
    TimIntImpl(const TimIntImpl& old);

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    virtual void set_element_time_parameter() const = 0;

    //! set time for evaluation of Neumann boundary conditions
    virtual void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) = 0;

    //! Set general element parameters
    void set_element_general_parameters() const;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/


    /*--- calculate and update -----------------------------------------------*/

    //! Apply Dirichlet boundary conditions on provided state vector
    void apply_dirichlet_bc(const double time,                //!< evaluation time
        std::shared_ptr<Core::LinAlg::Vector<double>> prenp,  //!< pressure (may be = null)
        std::shared_ptr<Core::LinAlg::Vector<double>>
            predt  //!< first time derivative (may be = null)
    );

    //! potential residual scaling and potential addition of Neumann terms
    void scaling_and_neumann();

    //! add actual Neumann loads multiplied with time factor to the residual
    virtual void add_neumann_to_residual() = 0;

    //! Apply Neumann boundary conditions
    void apply_neumann_bc(Core::LinAlg::Vector<double>& neumann_loads  //!< Neumann loads
    );

    //! call elements to calculate system matrix and rhs and assemble
    virtual void assemble_mat_and_rhs();

    //! return the right time-scaling-factor for the true residual
    virtual double residual_scaling() const = 0;

    //! penalty term to ensure positive pressures (cavitation)
    virtual void add_cavitation_penalty();

    //! contains the nonlinear iteration loop
    virtual void nonlinear_solve();

    //! check convergence (or divergence) of nonlinear iteration
    bool abort_nonlin_iter(const int itnum,  //!< current value of iteration step counter
        const int itemax,                    //!< maximum number of iteration steps
        const double ittol,                  //!< relative tolerance for increments
        const double abstolres,              //!< absolute tolerance for the residual norm
        double& actresidual                  //!< return value of the current residual
    );

    //! Calculate problem specific norm
    virtual void calc_problem_specific_norm(
        double& preresnorm, double& incprenorm_L2, double& prenorm_L2, double& preresnorminf);

    /*--- query and output ---------------------------------------------------*/

    //! is output needed for the current time step?
    bool do_output() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    //! write state vectors prenp to BINIO
    virtual void output_state();

    //! write state vectors prenp to Gmsh postprocessing files
    void output_to_gmsh(const int step, const double time) const;

    //! print header of convergence table to screen
    virtual void print_convergence_header();

    //! print first line of convergence table to screen
    virtual void print_convergence_values_first_iter(
        const int& itnum,            //!< current Newton-Raphson iteration step
        const int& itemax,           //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,         //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm,    //!< L2 norm of pressure residual
        const double& preresnorminf  //!< infinity norm of pressure residual
    );

    //! print current line of convergence table to screen
    virtual void print_convergence_values(
        const int& itnum,             //!< current Newton-Raphson iteration step
        const int& itemax,            //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,          //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm,     //!< L2 norm of pressure residual
        const double& incprenorm_L2,  //!< L2 norm of pressure increment
        const double& prenorm_L2,     //!< L2 norm of pressure state vector
        const double& preresnorminf   //!< infinity norm of pressure residual
    );

    //! print finish line of convergence table to screen
    virtual void print_convergence_finish_line();

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! increment time and step value
    void increment_time_and_step();

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! solver
    std::shared_ptr<Core::LinAlg::Solver> solver_;

    //! parameter list
    const std::shared_ptr<Teuchos::ParameterList> params_;

    //! processor id
    int myrank_;

    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for Eulerian or ALE formulation of equation(s)
    bool isale_;

    //! incremental or linear full solving? rename -> is_incremental_
    bool incremental_;

    //! flag for Modified Reynolds Equation
    bool modified_reynolds_;

    //! flag for adding squeeze term to Reynolds Equ.
    bool addsqz_;

    //! flag for pure lubrication problem
    bool purelub_;


    /*--- query and output ---------------------------------------------------*/

    //! flag for printing out mean values of pressures
    const bool outmean_;

    //! boolean to write Gmsh postprocessing files (input parameter)
    const bool outputgmsh_;

    //! boolean to write state vector to matlab file (input parameter)
    const bool output_state_matlab_;

    /*========================================================================*/
    //! @name Time, time-step, and iteration variables
    /*========================================================================*/

    //! actual time
    double time_;

    //! maximum simulation time
    double maxtime_;

    //! actual step number
    int step_;

    //! maximum number of steps ? name maxtime vs. stepmax
    int stepmax_;

    //! time step size
    double dta_;

    //! time measurement element
    double dtele_;

    //! time measurement solve
    double dtsolve_;

    //! number of newton iterations in actual timestep
    int iternum_;

    /*========================================================================*/
    //! @name pressure degrees of freedom variables
    /*========================================================================*/

    //! number of space dimensions
    int nsd_;

    //! pressure at time n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> prenp_;

    /*========================================================================*/
    //! @name velocity, pressure, and related
    /*========================================================================*/

    //! number of dofset associated with displacement dofs
    int nds_disp_;

    /*========================================================================*/
    //! @name Galerkin discretization, boundary conditions, and related
    /*========================================================================*/

    //! the lubrication discretization
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! the discretization writer
    std::shared_ptr<Core::IO::DiscretizationWriter> output_;

    //! system matrix (either sparse matrix or block sparse matrix)
    std::shared_ptr<Core::LinAlg::SparseOperator> sysmat_;

    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    std::shared_ptr<Core::LinAlg::Vector<double>> zeros_;

    //! maps for extracting Dirichlet and free DOF sets
    std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps_;

    //! the vector containing body and surface forces
    std::shared_ptr<Core::LinAlg::Vector<double>> neumann_loads_;

    //! residual vector
    std::shared_ptr<Core::LinAlg::Vector<double>> residual_;

    //! true (rescaled) residual vector without zeros at Dirichlet conditions
    std::shared_ptr<Core::LinAlg::Vector<double>> trueresidual_;

    //! nonlinear iteration increment vector
    std::shared_ptr<Core::LinAlg::Vector<double>> increment_;

    std::shared_ptr<Core::LinAlg::Vector<double>> prei_;  //!< residual pressures
                                                          //!< \f$\Delta{p}^{<k>}_{n+1}\f$

    //! Dirchlet toggle vector for unprojectable nodes (i.e. infinite gap)
    std::shared_ptr<Core::LinAlg::Vector<double>> inf_gap_toggle_lub_;

    /*========================================================================*/
    //! @name not classified variables - to be kept clean!!!
    /*========================================================================*/

    //! write results every upres_ steps ? writesolutionevery_
    int upres_;

    //! write restart data every uprestart_ steps ? writesolutioneveryrestart_
    int uprestart_;

    //! Surface roughness standard deviation used in Modified Reynolds Equation
    double roughness_deviation_;

    /*========================================================================*/


  };  // class TimIntImpl
}  // namespace Lubrication


FOUR_C_NAMESPACE_CLOSE

#endif
