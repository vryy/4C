/*----------------------------------------------------------------------*/
/*! \file

\brief Base class functions for time integration of electromagnetics

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef ELEMAG_TIMINT_H
#define ELEMAG_TIMINT_H

/*----------------------------------------------------------------------*
 | headers                                             gravemeier 06/17 |
 *----------------------------------------------------------------------*/
#include <fstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_inpar/inpar_elemag.H"

/*----------------------------------------------------------------------*
 | forward declarations                                gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseOperator;
  class MapExtractor;
  class Solver;
}  // namespace LINALG
namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class Node;
  class ResultTest;
}  // namespace DRT
namespace IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*
 | general time integration for electromagnetics       gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace ELEMAG
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
    ElemagTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// Virtual destructor.
    virtual ~ElemagTimeInt();

    /// Initialization routine.
    virtual void Init();

    /*!
    \brief Prints information about the discretization and time integration to screen.
    */
    virtual void PrintInformationToScreen();

    /*!
    \brief print the name of the scheme as std::string.
    */
    virtual std::string Name()
    {
      switch (elemagdyna_)
      {
        case INPAR::ELEMAG::DynamicType::elemag_bdf1:
          return "BDF1";
        // break;
        case INPAR::ELEMAG::DynamicType::elemag_bdf2:
          return "BDF2";
        case INPAR::ELEMAG::DynamicType::elemag_bdf4:
          return "BDF4";
          // break;
        default:
          return "General time integration class. No scheme implemented.";
      }
    };

    /*!
    \brief Print system matrix.
    */
    virtual void PrintSysmat_() { std::cout << sysmat_ << std::endl; };

    /*!
    \brief Iterates in time untill either the max number of steps or the final time has been
    reached.
    */
    virtual void Integrate();
    // void Integrate();

    /// Elements initialization
    virtual void ElementsInit();

    /*!
    \brief Set initial field by given function or null function.
    */
    virtual void SetInitialField(
        const INPAR::ELEMAG::InitialField init, Teuchos::ParameterList& start_params);

    /*!
    \brief Compare the numerical solution to the analytical one.
    */
    virtual Teuchos::RCP<Epetra_SerialDenseVector> ComputeError();

    /*!
    \brief Print the computed errors to screen.
    */
    virtual void PrintErrors(Teuchos::RCP<Epetra_SerialDenseVector>& errors);
    /*!
    \brief ProjectfieldTest is used for debugging purposes.

    This function projects a defined function to the interior fields and can be used to check if the
    global solver works. It is the dual of ProjectFieldTestTrace.
    */
    virtual void ProjectFieldTest(int startfuncno);

    /*!
    \brief ProjectfieldTestTrace is used for debugging purposes.

    This function projects a defined function to the global trace and can be used to check if the
    local solver works. It is the dual of ProjectFieldTest.
    */
    virtual void ProjectFieldTestTrace(int startfuncno);

    /*!
    \brief Initialize the algorithm using BDF1 as first step
    */
    void InitializeAlgorithm();

    /*!
    \brief Call elements to calculate system matrix and RHS.
    */
    virtual void AssembleMatAndRHS();
    // void AssembleMatAndRHS();

    /*!
    \brief Updates the local solution and updates the RHS wioth the new values.

    Once the solution of the global system has been carried on it is necessary to solve for the
    local problems and update the residual.
    */
    virtual void UpdateInteriorVariablesAndAssembleRHS();
    // void UpdateInteriorVariablesAndAssembleRHS();

    /*!
    \brief Apply Dirichlet boudnary conditions to system.
    */
    void ApplyDirichletToSystem();

    /*!
    \brief Compute Silver-Mueller boundary conditions.

    The functions ask the elements to compute the element matrices and vectors and assembles.

    \param[in]  resonly Boolean indicating if it is only necessary to recompute the RHS (true) or if
    it necessary to compute the matrices AND the RHS (false).
    */
    void ComputeSilverMueller(bool resonly);

    /*!
    \brief Solve the system for trace and interior variables.
    */
    virtual void Solve();
    // void Solve();

    /*!
    \brief Output to screen.
    */
    virtual void OutputToScreen();
    // void OutputToScreen();

    /*!
    \brief Output routine

    This routine outputs the computed global and local solutions to the output ".dat" file.
    */
    virtual void Output();
    // void Output();

    /*!
    \brief Write restart.

    This function allow the restart of a simulation from the result of the current
    simulation. It is needed to write the restart values that will be read by the #ReadRestart
    function.
    */
    virtual void WriteRestart();
    // void WriteRestart();

    /*!
    \brief Read restart.

    This function allow the restart of a simulation from the result of a previous
    simulation. It is needed to read the restart values that are written by the #WriteRestart
    function.
    */
    virtual void ReadRestart(int step);
    // void ReadRestart(int step);

    /*!
    \brief Obtain system matrix
    */
    Epetra_CrsMatrix SystemMatrix()
    {
      return (*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_, true)->EpetraMatrix());
    };

    /*!
    \brief print the system matrix giving the output filename

    */
    void SpySysmat(std::string filename, int precision = 20)
    {
      std::ofstream out;
      out.open(filename);
      SpySysmat(out << std::setprecision(precision));
      out.close();
    };
    /*!
    \brief print the system matrix to a given output stream

    This function still ha to be implemented.
    This function is supposed to work as a matlab or numpy spy().
    */
    void SpySysmat(std::ostream& out = std::cout);

    /*!
    \brief increment time and step value
    */
    virtual void IncrementTimeAndStep()
    {
      step_ += 1;
      time_ += dtp_;
    }

    double Time() { return time_; }
    int Step() { return step_; }
    int UpRes() { return upres_; }
    double TimeStep() { return dtp_; }

    /*!
    \brief returns pointer to the discretization
    */
    Teuchos::RCP<DRT::Discretization> Discretization();

    /*!
    \brief create result test
    */
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

   protected:
    /// discretization, solver, parameter list and output
    Teuchos::RCP<DRT::DiscretizationHDG> discret_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    INPAR::ELEMAG::DynamicType elemagdyna_;  /// time integration scheme

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
    int errfunct_;  /// Function number for error calculation

    int sourcefuncno_;  /// Source function number

    /// system matrix
    Teuchos::RCP<LINALG::SparseOperator> sysmat_;

    /// residual vector
    Teuchos::RCP<Epetra_Vector> residual_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

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

}  // namespace ELEMAG

#endif /* ELEMAG_TIMINT_H */
