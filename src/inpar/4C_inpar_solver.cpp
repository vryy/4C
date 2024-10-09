/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for linear solvers

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_solver.hpp"

#include "4C_linear_solver_method.hpp"
#include "4C_utils_parameter_list.hpp"

#include <BelosTypes.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar::SOLVER
{
  void set_valid_solver_parameters(Teuchos::ParameterList& list)
  {
    // Solver options
    {
      Teuchos::setStringToIntegralParameter<Core::LinearSolver::SolverType>("SOLVER", "undefined",
          "The solver to attack the system of linear equations arising of FE approach with.",
          Teuchos::tuple<std::string>("UMFPACK", "Superlu", "Belos", "undefined"),
          Teuchos::tuple<Core::LinearSolver::SolverType>(Core::LinearSolver::SolverType::umfpack,
              Core::LinearSolver::SolverType::superlu, Core::LinearSolver::SolverType::belos,
              Core::LinearSolver::SolverType::undefined),
          &list);
    }

    // Iterative solver options
    {
      Teuchos::setStringToIntegralParameter<Core::LinearSolver::IterativeSolverType>("AZSOLVE",
          "GMRES", "Type of linear solver algorithm to use.",
          Teuchos::tuple<std::string>("CG", "GMRES", "BiCGSTAB"),
          Teuchos::tuple<Core::LinearSolver::IterativeSolverType>(
              Core::LinearSolver::IterativeSolverType::cg,
              Core::LinearSolver::IterativeSolverType::gmres,
              Core::LinearSolver::IterativeSolverType::bicgstab),
          &list);
    }

    // Preconditioner options
    {
      Teuchos::setStringToIntegralParameter<Core::LinearSolver::PreconditionerType>("AZPREC", "ILU",
          "Type of internal preconditioner to use.\n"
          "Note! this preconditioner will only be used if the input operator\n"
          "supports the Epetra_RowMatrix interface and the client does not pass\n"
          "in an external preconditioner!",
          Teuchos::tuple<std::string>(
              "ILU", "ML", "MLFLUID2", "MueLu", "MueLu_contactSP", "AMGnxn", "CheapSIMPLE", "Teko"),
          Teuchos::tuple<Core::LinearSolver::PreconditionerType>(
              Core::LinearSolver::PreconditionerType::ilu,
              Core::LinearSolver::PreconditionerType::multigrid_ml,
              Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2,
              Core::LinearSolver::PreconditionerType::multigrid_muelu,
              Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp,
              Core::LinearSolver::PreconditionerType::multigrid_nxn,
              Core::LinearSolver::PreconditionerType::cheap_simple,
              Core::LinearSolver::PreconditionerType::block_teko),
          &list);
    }

    // Ifpack options
    {
      Core::UTILS::int_parameter("IFPACKOVERLAP", 0,
          "The amount of overlap used for the ifpack \"ilu\" preconditioner.", &list);

      Core::UTILS::int_parameter("IFPACKGFILL", 0,
          "The amount of fill allowed for an internal \"ilu\" preconditioner.", &list);

      std::vector<std::string> ifpack_combine_valid_input = {"Add", "Insert", "Zero"};
      Core::UTILS::string_parameter("IFPACKCOMBINE", "Add",
          "Combine mode for Ifpack Additive Schwarz", &list, ifpack_combine_valid_input);
    }

    // Iterative solver options
    {
      Core::UTILS::int_parameter("AZITER", 1000,
          "The maximum number of iterations the underlying iterative solver is allowed to "
          "perform",
          &list);

      Core::UTILS::double_parameter("AZTOL", 1e-8,
          "The level the residual norms must reach to decide about successful convergence", &list);

      Teuchos::setStringToIntegralParameter<Belos::ScaleType>("AZCONV", "AZ_r0",
          "The implicit residual norm scaling type to use for terminating the iterative solver.",
          Teuchos::tuple<std::string>("AZ_r0", "AZ_noscaled"),
          Teuchos::tuple<Belos::ScaleType>(Belos::ScaleType::NormOfInitRes, Belos::ScaleType::None),
          &list);

      Core::UTILS::int_parameter("AZOUTPUT", 0,
          "The number of iterations between each output of the solver's progress is written to "
          "screen",
          &list);
      Core::UTILS::int_parameter(
          "AZREUSE", 0, "The number specifying how often to recompute some preconditioners", &list);

      Core::UTILS::int_parameter("AZSUB", 50,
          "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
          "a restart is performed.",
          &list);

      Core::UTILS::int_parameter("AZGRAPH", 0, "unused", &list);
      Core::UTILS::int_parameter("AZBDIAG", 0, "unused", &list);
      Core::UTILS::double_parameter("AZOMEGA", 0.0, "unused", &list);

      Core::UTILS::string_parameter(
          "SOLVER_XML_FILE", "none", "xml file defining any linear solver", &list);
    }

    // ML options
    {
      Core::UTILS::int_parameter("ML_PRINT", 0, "ML print-out level (0-10)", &list);
      Core::UTILS::int_parameter(
          "ML_MAXCOARSESIZE", 5000, "ML stop coarsening when coarse ndof smaller then this", &list);
      Core::UTILS::int_parameter("ML_MAXLEVEL", 5, "ML max number of levels", &list);
      Core::UTILS::int_parameter("ML_AGG_SIZE", 27,
          "objective size of an aggregate with METIS/VBMETIS, 2D: 9, 3D: 27", &list);

      Core::UTILS::double_parameter("ML_DAMPFINE", 1., "damping fine grid", &list);
      Core::UTILS::double_parameter("ML_DAMPMED", 1., "damping med grids", &list);
      Core::UTILS::double_parameter("ML_DAMPCOARSE", 1., "damping coarse grid", &list);
      Core::UTILS::double_parameter("ML_PROLONG_SMO", 0.,
          "damping factor for prolongator smoother (usually 1.33 or 0.0)", &list);
      Core::UTILS::double_parameter(
          "ML_PROLONG_THRES", 0., "threshold for prolongator smoother/aggregation", &list);

      Core::UTILS::string_parameter("ML_SMOTIMES", "1 1 1 1 1",
          "no. smoothing steps or polynomial order on each level (at least ML_MAXLEVEL numbers)",
          &list);

      std::vector<std::string> ml_coarsen_valid_input = {"UC", "METIS", "VBMETIS", "MIS"};
      Core::UTILS::string_parameter("ML_COARSEN", "UC", "", &list, ml_coarsen_valid_input);

      Core::UTILS::bool_parameter("ML_REBALANCE", "Yes",
          "Performe ML-internal rebalancing of coarse level operators.", &list);

      std::vector<std::string> ml_smootherfine_valid_input = {"SGS", "Jacobi", "Chebychev", "MLS",
          "ILU", "KLU", "Superlu", "GS", "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD",
          "Uzawa"};
      Core::UTILS::string_parameter(
          "ML_SMOOTHERFINE", "ILU", "", &list, ml_smootherfine_valid_input);

      std::vector<std::string> ml_smoothermed_valid_input = {"SGS", "Jacobi", "Chebychev", "MLS",
          "ILU", "KLU", "Superlu", "GS", "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD",
          "Uzawa"};
      Core::UTILS::string_parameter("ML_SMOOTHERMED", "ILU", "", &list, ml_smoothermed_valid_input);

      std::vector<std::string> ml_smoothercoarse_valid_input = {"SGS", "Jacobi", "Chebychev", "MLS",
          "ILU", "KLU", "Superlu", "GS", "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD",
          "Uzawa"};
      Core::UTILS::string_parameter(
          "ML_SMOOTHERCOARSE", "Umfpack", "", &list, ml_smoothercoarse_valid_input);


      Core::UTILS::int_parameter("SUB_SOLVER1", -1,
          "sub solver/smoother block number (SIMPLE/C: used for prediction of primary variable "
          "on "
          "all levels, BS: used for fine and intermedium BraessSarazin (BS) level smoother)",
          &list);
      Core::UTILS::int_parameter("SUB_SOLVER2", -1,
          "sub solver/smoother block number (SIMPLE/C: used for SchurComplement eq. on all "
          "levels, "
          "BS: used for coarse BraessSarazin (BS) level smoother)",
          &list);
    }

    // MueLu options
    {
      Core::UTILS::string_parameter(
          "MUELU_XML_FILE", "none", "xml file defining any MueLu preconditioner", &list);
    }

    // Teko options
    {
      Core::UTILS::string_parameter(
          "TEKO_XML_FILE", "none", "xml file defining any Teko preconditioner", &list);
    }

    // user-given name of solver block (just for beauty)
    Core::UTILS::string_parameter("NAME", "No_name", "User specified name for solver block", &list);

    // damping parameter for SIMPLE
    Core::UTILS::double_parameter(
        "SIMPLE_DAMPING", 1., "damping parameter for SIMPLE preconditioner", &list);

    // Parameters for AMGnxn Preconditioner
    {
      Core::UTILS::string_parameter("AMGNXN_TYPE", "AMG(BGS)",
          "Name of the pre-built preconditioner to be used. If set to\"XML\" the preconditioner "
          "is defined using a xml file",
          &list);
      Core::UTILS::string_parameter(
          "AMGNXN_XML_FILE", "none", "xml file defining the AMGnxn preconditioner", &list);
    }
  }


  void set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
  {
    // set valid parameters for solver blocks

    // Note: the maximum number of solver blocks is hardwired here. If you change this,
    // don't forget to edit the corresponding parts in globalproblems.cpp, too.
    for (int i = 1; i < 10; i++)
    {
      std::stringstream ss;
      ss << "SOLVER " << i;
      std::stringstream ss_description;
      ss_description << "solver parameters for solver block " << i;
      Teuchos::ParameterList& solverlist = list->sublist(ss.str(), false, ss_description.str());
      set_valid_solver_parameters(solverlist);
    }
  }

}  // namespace Inpar::SOLVER

FOUR_C_NAMESPACE_CLOSE
