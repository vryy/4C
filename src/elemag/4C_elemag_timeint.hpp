/*----------------------------------------------------------------------*/
/*! \file

\brief Base class functions for time integration of electromagnetics

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_TIMEINT_HPP
#define FOUR_C_ELEMAG_TIMEINT_HPP

/*----------------------------------------------------------------------*
 | headers                                             gravemeier 06/17 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_elemag.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace Core::LinAlg
{
  class SparseOperator;
  class EquilibrationSparse;
  enum class EquilibrationMethod;
  class MapExtractor;
  class Solver;
}  // namespace Core::LinAlg
namespace Core::FE
{
  class Discretization;
  class DiscretizationHDG;
}  // namespace Core::FE

namespace Core::IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*
 | general time integration for electromagnetics       gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace EleMag
{
  /*!
  \brief Base class functions for time integration of electromagnetics

  This has to be the base class that is inherited by all the specific time integration rules but for
  now it includes all the time integration. This will be fixed in the following iteration of the
  code itself.

  \author Berardocco
  */
  class ElemagTimeInt
  {
   public:
    /// Constructor.
    ElemagTimeInt(const Teuchos::RCP<Core::FE::DiscretizationHDG>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output);

    /// Virtual destructor.
    virtual ~ElemagTimeInt() = default;

    /// Initialization routine.
    virtual void init();

    /*!
    \brief Prints information about the discretization and time integration to screen.
    */
    virtual void print_information_to_screen();

    /*!
    \brief print the name of the scheme as std::string.
    */
    virtual std::string name()
    {
      switch (elemagdyna_)
      {
        case Inpar::EleMag::DynamicType::elemag_bdf1:
          return "BDF1";
        // break;
        case Inpar::EleMag::DynamicType::elemag_bdf2:
          return "BDF2";
        case Inpar::EleMag::DynamicType::elemag_bdf4:
          return "BDF4";
          // break;
        default:
          return "General time integration class. No scheme implemented.";
      }
    };

    /*!
    \brief Print system matrix.
    */
    virtual void print_sysmat() { std::cout << sysmat_ << std::endl; };

    /*!
    \brief Iterates in time untill either the max number of steps or the final time has been
    reached.
    */
    virtual void integrate();
    // void Integrate();

    /// Elements initialization
    virtual void elements_init();

    /*!
    \brief Set initial field by given function or null function.
    */
    virtual void set_initial_field(const Inpar::EleMag::InitialField init, int startfuncno);

    /*!
    \brief Import initial electric field from scatra solution
    */
    void set_initial_electric_field(
        Teuchos::RCP<Epetra_Vector> phi, Teuchos::RCP<Core::FE::Discretization>& scatradis);

    /*!
    \brief Compare the numerical solution to the analytical one.
    */
    virtual Teuchos::RCP<Core::LinAlg::SerialDenseVector> compute_error();

    /*!
    \brief Print the computed errors to screen.
    */
    virtual void print_errors(Teuchos::RCP<Core::LinAlg::SerialDenseVector>& errors);
    /*!
    \brief ProjectfieldTest is used for debugging purposes.

    This function projects a defined function to the interior fields and can be used to check if the
    global solver works. It is the dual of project_field_test_trace.
    */
    virtual void project_field_test(int startfuncno);

    /*!
    \brief ProjectfieldTestTrace is used for debugging purposes.

    This function projects a defined function to the global trace and can be used to check if the
    local solver works. It is the dual of ProjectFieldTest.
    */
    virtual void project_field_test_trace(int startfuncno);

    /*!
    \brief Initialize the algorithm using BDF1 as first step
    */
    void initialize_algorithm();

    /*!
    \brief Call elements to calculate system matrix and RHS.
    */
    virtual void assemble_mat_and_rhs();
    // void assemble_mat_and_rhs();

    /*!
    \brief Updates the local solution and updates the RHS wioth the new values.

    Once the solution of the global system has been carried on it is necessary to solve for the
    local problems and update the residual.
    */
    virtual void update_interior_variables_and_assemble_rhs();
    // void update_interior_variables_and_assemble_rhs();

    /*!
    \brief Apply Dirichlet boudnary conditions to system.

    \param[in]  resonly Boolean indicating if it is only necessary to recompute the RHS (true) or if
    it necessary to compute the matrices AND the RHS (false).
    */
    void apply_dirichlet_to_system(bool resonly);

    /*!
    \brief Compute Silver-Mueller boundary conditions.

    The functions ask the elements to compute the element matrices or vectors related to the
    Silver--Mueller (absorbing) boundary conditions and assembles.

    \param[in]  do_rhs Boolean indicating if we want to compute the rhs OR the matrix.

    \note We can not do both the matrix and the rhs as the underlying evaluate_condition does not
    replace the values in the rhs vector but assembles it. Therefore, by doing both the matrix and
    the vector we would end up summing more than one contribution at once
    */
    void compute_silver_mueller(bool do_rhs);

    /*!
    \brief Solve the system for trace and interior variables.
    */
    virtual void solve();
    // void Solve();

    /*!
    \brief Output to screen.
    */
    virtual void output_to_screen();
    // void OutputToScreen();

    /*!
    \brief Output routine

    This routine outputs the computed global and local solutions to the output ".dat" file.
    */
    virtual void output();
    // void output();

    /*!
    \brief Write restart.

    This function allow the restart of a simulation from the result of the current
    simulation. It is needed to write the restart values that will be read by the #read_restart
    function.
    */
    virtual void write_restart();
    // void write_restart();

    /*!
    \brief Read restart.

    This function allow the restart of a simulation from the result of a previous
    simulation. It is needed to read the restart values that are written by the #write_restart
    function.
    */
    virtual void read_restart(int step);
    // void read_restart(int step);

    /*!
    \brief Obtain system matrix
    */
    Epetra_CrsMatrix system_matrix()
    {
      return (
          *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_, true)->epetra_matrix());
    };

    /*!
    \brief print the system matrix giving the output filename

    */
    void spy_sysmat(std::string filename, int precision = 20)
    {
      std::ofstream out;
      out.open(filename);
      spy_sysmat(out << std::setprecision(precision));
      out.close();
    };
    /*!
    \brief print the system matrix to a given output stream

    This function still ha to be implemented.
    This function is supposed to work as a matlab or numpy spy().
    */
    void spy_sysmat(std::ostream& out = std::cout);

    /*!
    \brief increment time and step value
    */
    virtual void increment_time_and_step()
    {
      step_ += 1;
      time_ += dtp_;
    }

    double time() { return time_; }
    int step() { return step_; }
    int up_res() { return upres_; }
    double time_step() { return dtp_; }

    /*!
    \brief returns pointer to the discretization
    */
    Teuchos::RCP<Core::FE::Discretization> discretization();

    /*!
    \brief create result test
    */
    virtual Teuchos::RCP<Core::UTILS::ResultTest> create_field_test();

   protected:
    /// discretization, solver, parameter list and output
    Teuchos::RCP<Core::FE::DiscretizationHDG> discret_;
    Teuchos::RCP<Core::LinAlg::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<Core::IO::DiscretizationWriter> output_;

    Inpar::EleMag::DynamicType elemagdyna_;  /// time integration scheme

    int myrank_;  /// processor id

    double time_;  /// physical time
    int step_;     /// time step
    int restart_;  /// restart step

    double maxtime_;  /// maximum time
    int stepmax_;     /// maximum step

    int uprestart_;  /// write restart data every uprestart_ steps
    int upres_;      /// write output every upres_ steps

    int numdim_;        /// number of spatial dimensions
    double dtp_;        /// time step size
    const double tau_;  /// stabilization parameter

    double dtele_;    /// element evaluation time
    double dtsolve_;  /// solver time

    bool calcerr_;  /// flag for error calculation

    bool postprocess_;  /// flag for postprocessing

    int errfunct_;  /// Function number for error calculation

    int sourcefuncno_;  /// Source function number

    Core::LinAlg::EquilibrationMethod equilibration_method_;  /// equilibration method

    /// system matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat_;

    /// residual vector
    Teuchos::RCP<Epetra_Vector> residual_;

    /// all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<Core::LinAlg::EquilibrationSparse> equilibration_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;

    /// vector of zeros to be used for enforcing zero Dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;


    //  The fomulation relies on the variables:
    //  o   E (electri field)
    //  o   H (magnetic field field)
    //  o   \Lambda (hybrid varible)

    /// Trace vector to be solved at every iteration
    Teuchos::RCP<Epetra_Vector> trace_;

    /// Output vectors
    Teuchos::RCP<Epetra_MultiVector> electric;
    Teuchos::RCP<Epetra_MultiVector> electric_post;
    Teuchos::RCP<Epetra_MultiVector> magnetic;
    Teuchos::RCP<Epetra_MultiVector> trace;
    Teuchos::RCP<Epetra_Vector> conductivity;
    Teuchos::RCP<Epetra_Vector> permittivity;
    Teuchos::RCP<Epetra_Vector> permeability;
  };

}  // namespace EleMag

FOUR_C_NAMESPACE_CLOSE

#endif
