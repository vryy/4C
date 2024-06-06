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
  void SetValidSolverParameters(Teuchos::ParameterList& list)
  {
    using Teuchos::setStringToIntegralParameter;
    using Teuchos::tuple;

    // Solver options
    {
      setStringToIntegralParameter<Core::LinearSolver::SolverType>("SOLVER", "undefined",
          "The solver to attack the system of linear equations arising of FE approach with.",
          tuple<std::string>("UMFPACK", "Superlu", "Belos", "undefined"),
          tuple<Core::LinearSolver::SolverType>(Core::LinearSolver::SolverType::umfpack,
              Core::LinearSolver::SolverType::superlu, Core::LinearSolver::SolverType::belos,
              Core::LinearSolver::SolverType::undefined),
          &list);
    }

    // Iterative solver options
    {
      setStringToIntegralParameter<Core::LinearSolver::IterativeSolverType>("AZSOLVE", "GMRES",
          "Type of linear solver algorithm to use.", tuple<std::string>("CG", "GMRES", "BiCGSTAB"),
          tuple<Core::LinearSolver::IterativeSolverType>(
              Core::LinearSolver::IterativeSolverType::cg,
              Core::LinearSolver::IterativeSolverType::gmres,
              Core::LinearSolver::IterativeSolverType::bicgstab),
          &list);
    }

    // Preconditioner options
    {
      // this one is longer than 15 and the tuple<> function does not support this,
      // so build the Tuple class directly (which can be any size)
      Teuchos::Tuple<std::string, 14> name;
      Teuchos::Tuple<Core::LinearSolver::PreconditionerType, 14> number;

      name[0] = "ILU";
      number[0] = Core::LinearSolver::PreconditionerType::ilu;
      name[1] = "ICC";
      number[1] = Core::LinearSolver::PreconditionerType::icc;
      name[2] = "ML";
      number[2] = Core::LinearSolver::PreconditionerType::multigrid_ml;
      name[3] = "MLFLUID";
      number[3] = Core::LinearSolver::PreconditionerType::multigrid_ml_fluid;
      name[4] = "MLFLUID2";
      number[4] = Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2;
      name[5] = "MueLu";
      number[5] = Core::LinearSolver::PreconditionerType::multigrid_muelu;
      name[6] = "MueLu_fluid";
      number[6] = Core::LinearSolver::PreconditionerType::multigrid_muelu_fluid;
      name[7] = "MueLu_tsi";
      number[7] = Core::LinearSolver::PreconditionerType::multigrid_muelu_tsi;
      name[8] = "MueLu_contactSP";
      number[8] = Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp;
      name[9] = "MueLu_BeamSolid";
      number[9] = Core::LinearSolver::PreconditionerType::multigrid_muelu_beamsolid;
      name[10] = "MueLu_fsi";
      number[10] = Core::LinearSolver::PreconditionerType::multigrid_muelu_fsi;
      name[11] = "AMGnxn";
      number[11] = Core::LinearSolver::PreconditionerType::multigrid_nxn;
      name[12] = "BGS2x2";
      number[12] = Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2;
      name[13] = "CheapSIMPLE";
      number[13] = Core::LinearSolver::PreconditionerType::cheap_simple;

      setStringToIntegralParameter<Core::LinearSolver::PreconditionerType>("AZPREC", "ILU",
          "Type of internal preconditioner to use.\n"
          "Note! this preconditioner will only be used if the input operator\n"
          "supports the Epetra_RowMatrix interface and the client does not pass\n"
          "in an external preconditioner!",
          name, number, &list);
    }

    // Ifpack options
    {
      Core::UTILS::IntParameter("IFPACKOVERLAP", 0,
          "The amount of overlap used for the ifpack \"ilu\" preconditioner.", &list);

      Core::UTILS::IntParameter("IFPACKGFILL", 0,
          "The amount of fill allowed for an internal \"ilu\" preconditioner.", &list);

      setStringToIntegralParameter<int>("IFPACKCOMBINE", "Add",
          "Combine mode for Ifpack Additive Schwarz", tuple<std::string>("Add", "Insert", "Zero"),
          tuple<int>(0, 1, 2), &list);
    }

    // Iterative solver options
    {
      Core::UTILS::IntParameter("AZITER", 1000,
          "The maximum number of iterations the underlying iterative solver is allowed to perform",
          &list);

      Core::UTILS::DoubleParameter("AZTOL", 1e-8,
          "The level the residual norms must reach to decide about successful convergence", &list);

      setStringToIntegralParameter<Belos::ScaleType>("AZCONV", "AZ_r0",
          "The implicit residual norm scaling type to use for terminating the iterative solver.",
          tuple<std::string>("AZ_r0", "AZ_noscaled"),
          tuple<Belos::ScaleType>(Belos::ScaleType::NormOfInitRes, Belos::ScaleType::None), &list);

      Core::UTILS::IntParameter("AZOUTPUT", 0,
          "The number of iterations between each output of the solver's progress is written to "
          "screen",
          &list);
      Core::UTILS::IntParameter(
          "AZREUSE", 0, "The number specifying how often to recompute some preconditioners", &list);

      Core::UTILS::IntParameter("AZSUB", 50,
          "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
          "a restart is performed.",
          &list);

      Core::UTILS::IntParameter("AZGRAPH", 0, "unused", &list);
      Core::UTILS::IntParameter("AZBDIAG", 0, "unused", &list);
      Core::UTILS::DoubleParameter("AZOMEGA", 0.0, "unused", &list);
    }

    // ML options
    {
      Core::UTILS::IntParameter("ML_PRINT", 0, "ML print-out level (0-10)", &list);
      Core::UTILS::IntParameter(
          "ML_MAXCOARSESIZE", 5000, "ML stop coarsening when coarse ndof smaller then this", &list);
      Core::UTILS::IntParameter("ML_MAXLEVEL", 5, "ML max number of levels", &list);
      Core::UTILS::IntParameter("ML_AGG_SIZE", 27,
          "objective size of an aggregate with METIS/VBMETIS, 2D: 9, 3D: 27", &list);

      Core::UTILS::DoubleParameter("ML_DAMPFINE", 1., "damping fine grid", &list);
      Core::UTILS::DoubleParameter("ML_DAMPMED", 1., "damping med grids", &list);
      Core::UTILS::DoubleParameter("ML_DAMPCOARSE", 1., "damping coarse grid", &list);
      Core::UTILS::DoubleParameter("ML_PROLONG_SMO", 0.,
          "damping factor for prolongator smoother (usually 1.33 or 0.0)", &list);
      Core::UTILS::DoubleParameter(
          "ML_PROLONG_THRES", 0., "threshold for prolongator smoother/aggregation", &list);

      Core::UTILS::StringParameter("ML_SMOTIMES", "1 1 1 1 1",
          "no. smoothing steps or polynomial order on each level (at least ML_MAXLEVEL numbers)",
          &list);

      setStringToIntegralParameter<int>("ML_COARSEN", "UC", "",
          tuple<std::string>("UC", "METIS", "VBMETIS", "MIS"), tuple<int>(0, 1, 2, 3), &list);

      Core::UTILS::BoolParameter("ML_REBALANCE", "Yes",
          "Performe ML-internal rebalancing of coarse level operators.", &list);

      setStringToIntegralParameter<int>("ML_SMOOTHERFINE", "ILU", "",
          tuple<std::string>("SGS", "Jacobi", "Chebychev", "MLS", "ILU", "KLU", "Superlu", "GS",
              "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD", "Uzawa"),
          tuple<int>(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), &list);

      setStringToIntegralParameter<int>("ML_SMOOTHERMED", "ILU", "",
          tuple<std::string>("SGS", "Jacobi", "Chebychev", "MLS", "ILU", "KLU", "Superlu", "GS",
              "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD", "Uzawa"),
          tuple<int>(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), &list);

      setStringToIntegralParameter<int>("ML_SMOOTHERCOARSE", "Umfpack", "",
          tuple<std::string>("SGS", "Jacobi", "Chebychev", "MLS", "ILU", "KLU", "Superlu", "GS",
              "DGS", "Umfpack", "BS", "SIMPLE", "SIMPLEC", "IBD", "Uzawa"),
          tuple<int>(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), &list);

      Core::UTILS::IntParameter("SUB_SOLVER1", -1,
          "sub solver/smoother block number (SIMPLE/C: used for prediction of primary variable "
          "on "
          "all levels, BS: used for fine and intermedium BraessSarazin (BS) level smoother)",
          &list);
      Core::UTILS::IntParameter("SUB_SOLVER2", -1,
          "sub solver/smoother block number (SIMPLE/C: used for SchurComplement eq. on all "
          "levels, "
          "BS: used for coarse BraessSarazin (BS) level smoother)",
          &list);
    }

    // MueLu options
    {
      Core::UTILS::StringParameter(
          "MUELU_XML_FILE", "none", "xml file defining any MueLu preconditioner", &list);

      Core::UTILS::BoolParameter(
          "MUELU_XML_ENFORCE", "Yes", "option defining xml file usage", &list);
    }

    // BGS2x2 options
    {
      // switch order of blocks in BGS2x2 preconditioner
      setStringToIntegralParameter<int>("BGS2X2_FLIPORDER", "block0_block1_order",
          "BGS2x2 flip order parameter",
          tuple<std::string>("block0_block1_order", "block1_block0_order"), tuple<int>(0, 1),
          &list);

      // damping parameter for BGS2X2
      Core::UTILS::DoubleParameter(
          "BGS2X2_GLOBAL_DAMPING", 1., "damping parameter for BGS2X2 preconditioner", &list);
      Core::UTILS::DoubleParameter(
          "BGS2X2_BLOCK1_DAMPING", 1., "damping parameter for BGS2X2 preconditioner block1", &list);
      Core::UTILS::DoubleParameter(
          "BGS2X2_BLOCK2_DAMPING", 1., "damping parameter for BGS2X2 preconditioner block2", &list);
    }

    // user-given name of solver block (just for beauty)
    Core::UTILS::StringParameter("NAME", "No_name", "User specified name for solver block", &list);

    // damping parameter for SIMPLE
    Core::UTILS::DoubleParameter(
        "SIMPLE_DAMPING", 1., "damping parameter for SIMPLE preconditioner", &list);

    // Parameters for AMGnxn Preconditioner
    {
      Core::UTILS::StringParameter("AMGNXN_TYPE", "AMG(BGS)",
          "Name of the pre-built preconditioner to be used. If set to\"XML\" the preconditioner "
          "is defined using a xml file",
          &list);
      Core::UTILS::StringParameter(
          "AMGNXN_XML_FILE", "none", "xml file defining the AMGnxn preconditioner", &list);
    }
  }


  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
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
      SetValidSolverParameters(solverlist);
    }

    /*----------------------------------------------------------------------*/
    // UMFPACK solver section
    // some people just need a solver quickly. We provide a special parameter set
    // for UMFPACK that users can just use temporarily without introducing a
    // separate solver block.
    Teuchos::ParameterList& solver_u =
        list->sublist("UMFPACK SOLVER", false, "solver parameters for UMFPACK");
    SetValidSolverParameters(solver_u);
  }

}  // namespace Inpar::SOLVER

FOUR_C_NAMESPACE_CLOSE
