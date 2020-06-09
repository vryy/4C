/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for linear solvers

\maintainer Martin Kronbichler

\level 1
*/

/*----------------------------------------------------------------------*/

#include "inpar_solver.H"
#include "drt_validparameters.H"

#include <AztecOO.h>


namespace INPAR
{
  namespace SOLVER
  {
    void SetValidSolverParameters(Teuchos::ParameterList& list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;
      using namespace DRT::INPUT;

      Teuchos::Tuple<std::string, 12> solver_name;
      Teuchos::Tuple<int, 12> solver_number;

      solver_name[0] = "Amesos_KLU_sym";
      solver_number[0] = amesos_klu_sym;
      solver_name[1] = "Amesos_KLU_nonsym";
      solver_number[1] = amesos_klu_nonsym;
      solver_name[2] = "Superlu";
      solver_number[2] = superlu;
      solver_name[3] = "Aztec_MSR";
      solver_number[3] = aztec_msr;
      solver_name[4] = "LAPACK_sym";
      solver_number[4] = lapack_sym;
      solver_name[5] = "LAPACK_nonsym";
      solver_number[5] = lapack_nonsym;
      solver_name[6] = "UMFPACK";
      solver_number[6] = umfpack;
      solver_name[7] = "Belos";
      solver_number[7] = belos;
      solver_name[8] = "Stratimikos_Amesos";
      solver_number[8] = stratimikos_amesos;
      solver_name[9] = "Stratimikos_Aztec";
      solver_number[9] = stratimikos_aztec;
      solver_name[10] = "Stratimikos_Belos";
      solver_number[10] = stratimikos_belos;
      solver_name[11] = "undefined";
      solver_number[11] = undefined;


      setStringToIntegralParameter<int>("SOLVER", "undefined",
          "The solver to attack the system of linear equations arising of FE approach with.",
          solver_name, solver_number, &list);

      setStringToIntegralParameter<int>("AZSOLVE", "GMRES",
          "Type of linear solver algorithm to use.",
          tuple<std::string>("CG", "GMRES", "GMRESR", "CGS", "TFQMR", "BiCGSTAB", "LU", "FGMRES"),
          tuple<int>(azsolv_CG, azsolv_GMRES, azsolv_GMRESR, azsolv_CGS, azsolv_TFQMR,
              azsolv_BiCGSTAB, azsolv_LU, belos_FGMRES),
          &list);

      {
        // this one is longer than 15 and the tuple<> function does not support this,
        // so build the Tuple class directly (which can be any size)
        Teuchos::Tuple<std::string, 28> name;
        Teuchos::Tuple<int, 28> number;

        name[0] = "none";
        number[0] = azprec_none;
        name[1] = "ILU";
        number[1] = azprec_ILU;
        name[2] = "ILUT";
        number[2] = azprec_ILUT;
        name[3] = "Jacobi";
        number[3] = azprec_Jacobi;
        name[4] = "SymmGaussSeidel";
        number[4] = azprec_SymmGaussSeidel;
        name[5] = "Least_Squares";
        number[5] = azprec_Least_Squares;
        name[6] = "Neumann";
        number[6] = azprec_Neumann;
        name[7] = "ICC";
        number[7] = azprec_ICC;
        name[8] = "LU";
        number[8] = azprec_LU;
        name[9] = "RILU";
        number[9] = azprec_RILU;
        name[10] = "ML";
        number[10] = azprec_ML;
        name[11] = "MLFLUID";
        number[11] = azprec_MLfluid;
        name[12] = "MLFLUID2";
        number[12] = azprec_MLfluid2;
        name[13] = "MLAPI";
        number[13] = azprec_MLAPI;
        name[14] = "GaussSeidel";
        number[14] = azprec_GaussSeidel;
        name[15] = "DownwindGaussSeidel";
        number[15] = azprec_DownwindGaussSeidel;
        name[16] = "BGS2x2";
        number[16] = azprec_BGS2x2;
        name[17] = "BGSnxn";
        number[17] = azprec_BGSnxn;
        name[18] = "TekoSIMPLE";
        number[18] = azprec_TekoSIMPLE;
        name[19] = "CheapSIMPLE";
        number[19] = azprec_CheapSIMPLE;
        name[20] = "MueLu_sym";
        number[20] = azprec_MueLuAMG_sym;
        name[21] = "MueLu_nonsym";
        number[21] = azprec_MueLuAMG_nonsym;
        name[22] = "MueLu_contact";
        number[22] = azprec_MueLuAMG_contact;
        name[23] = "MueLu_contact2";
        number[23] = azprec_MueLuAMG_contact2;
        name[24] = "MueLu_contactSP";
        number[24] = azprec_MueLuAMG_contactSP;
        name[25] = "MueLu_contactPenalty";
        number[25] = azprec_MueLuAMG_contactPen;
        name[26] = "AMGnxn";
        number[26] = azprec_AMGnxn;
        name[27] = "Chebyshev";
        number[27] = azprec_Chebyshev;

        setStringToIntegralParameter<int>("AZPREC", "ILU",
            "Type of internal preconditioner to use.\n"
            "Note! this preconditioner will only be used if the input operator\n"
            "supports the Epetra_RowMatrix interface and the client does not pass\n"
            "in an external preconditioner!",
            name, number, &list);
      }

      IntParameter("IFPACKOVERLAP", 0,
          "The amount of overlap used for the ifpack \"ilu\" and \"ilut\" preconditioners.", &list);

      IntParameter("IFPACKGFILL", 0,
          "This parameter has multiple meanings:\n1. The amount of fill-in into the graph allowed "
          "for the internal \"ilu\" preconditioner\n2. Number of sweeos for IFPACK-based point "
          "relaxation schemes.",
          &list);

      DoubleParameter("IFPACKFILL", 1.0,
          "The amount of fill allowed for an internal \"ilut\" preconditioner.", &list);

      setStringToIntegralParameter<int>("IFPACKCOMBINE", "Add",
          "Combine mode for Ifpack Additive Schwarz", tuple<std::string>("Add", "Insert", "Zero"),
          tuple<int>(0, 1, 2), &list);

      DoubleParameter("AZDROP", 0.0,
          "The tolerance below which an entry from the factors of an internal \"ilut\"\n"
          "preconditioner will be dropped.",
          &list);

      IntParameter("AZPOLY", 3,
          "The order for of the polynomials used for the \"Polynomial\" and\n"
          "\"Least-squares Polynomial\" internal preconditioners.",
          &list);

      IntParameter("AZSUB", 50,
          "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
          "a restart is performed.",
          &list);
      setStringToIntegralParameter<int>("AZCONV", "AZ_r0",  // Same as "rhs" when x=0
          "The convergence test to use for terminating the iterative solver.",
          tuple<std::string>("AZ_r0", "AZ_rhs", "AZ_Anorm", "AZ_noscaled", "AZ_sol", "AZ_weighted",
              "AZ_expected_values", "AZTECOO_conv_test", "AZ_inf_noscaled"),
          tuple<int>(AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol, AZ_weighted, AZ_expected_values,
              AZTECOO_conv_test, AZ_inf_noscaled),
          &list);

      IntParameter("AZOUTPUT", 0,  // By default, no output from Aztec!
          "The number of iterations between each output of the solver's progress.", &list);

      IntParameter("AZREUSE", 0, "how often to recompute some preconditioners", &list);
      IntParameter("AZITER", 1000, "max iterations", &list);
      IntParameter("AZGRAPH", 0, "unused", &list);
      IntParameter("AZBDIAG", 0, "", &list);

      DoubleParameter("AZTOL", 1e-8, "tolerance in (un)scaled residual", &list);
      DoubleParameter("AZOMEGA", 0.0, "damping for GaussSeidel and jacobi type methods", &list);
      DoubleParameter("DWINDTAU", 1.5, "threshold tau for downwinding", &list);

      setStringToIntegralParameter<int>("AZSCAL", "none", "scaling of the system",
          tuple<std::string>("none", "sym", "infnorm"), tuple<int>(0, 1, 2), &list);

      // parameters of ML preconditioner

      IntParameter("ML_PRINT", 0, "ML print-out level (0-10)", &list);
      IntParameter(
          "ML_MAXCOARSESIZE", 5000, "ML stop coarsening when coarse ndof smaller then this", &list);
      IntParameter("ML_MAXLEVEL", 5, "ML max number of levels", &list);
      IntParameter("ML_AGG_SIZE", 27,
          "objective size of an aggregate with METIS/VBMETIS, 2D: 9, 3D: 27", &list);

      DoubleParameter("ML_DAMPFINE", 1., "damping fine grid", &list);
      DoubleParameter("ML_DAMPMED", 1., "damping med grids", &list);
      DoubleParameter("ML_DAMPCOARSE", 1., "damping coarse grid", &list);
      DoubleParameter("ML_PROLONG_SMO", 0.,
          "damping factor for prolongator smoother (usually 1.33 or 0.0)", &list);
      DoubleParameter(
          "ML_PROLONG_THRES", 0., "threshold for prolongator smoother/aggregation", &list);

      StringParameter("ML_SMOTIMES", "1 1 1 1 1",
          "no. smoothing steps or polynomial order on each level (at least ML_MAXLEVEL numbers)",
          &list);

      setStringToIntegralParameter<int>("ML_COARSEN", "UC", "",
          tuple<std::string>("UC", "METIS", "VBMETIS", "MIS"), tuple<int>(0, 1, 2, 3), &list);

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

      // TODO remove this
      setStringToIntegralParameter<int>("MueLu_INITSMOOTHER", "SGS", "",
          tuple<std::string>("SGS", "Jacobi", "Chebychev", "ILU", "GS"), tuple<int>(0, 1, 2, 4, 7),
          &list);

      IntParameter("SUB_SOLVER1", -1,
          "sub solver/smoother block number (SIMPLE/C: used for prediction of primary variable on "
          "all levels, BS: used for fine and intermedium BraessSarazin (BS) level smoother)",
          &list);
      IntParameter("SUB_SOLVER2", -1,
          "sub solver/smoother block number (SIMPLE/C: used for SchurComplement eq. on all levels, "
          "BS: used for coarse BraessSarazin (BS) level smoother)",
          &list);

      // TODO remove this
      IntParameter("MueLu_INITSMOO_SWEEPS", 1,
          "number of sweeps for adaptive SA smoother (initialization phase). For Chebyshev it is "
          "used as polynomial degree",
          &list);
      DoubleParameter("MueLu_INITSMOO_DAMPING", 1.,
          "damping parameter for adaptive SA smoother (initialization phase). For Chebyshev it is "
          "used as alpha parameter",
          &list);

      IntParameter("MueLu_MIN_AGG_SIZE", 6,
          "Minimal objective size of an aggregate (to influence the coarsening rate)", &list);

      setStringToIntegralParameter<int>("MueLu_REBALANCE", "No",
          "activate rebalancing using Zoltan/Isorropia",
          tuple<std::string>("NO", "No", "no", "YES", "Yes", "yes"), tuple<int>(0, 1, 2, 3, 4, 5),
          &list);
      DoubleParameter("MueLu_REBALANCE_NONZEROIMBALANCE", 1.2,
          "maximum allowed nonzero imbalance factor", &list);
      IntParameter("MueLu_REBALANCE_MINROWS", 1000,
          "minimum numbers of rows per processor before rebalancing is necessary", &list);

      {
        // MueLu_Reuse
        Teuchos::Tuple<std::string, 3> name;
        Teuchos::Tuple<int, 3> number;

        name[0] = "nothing";
        number[0] = Reuse_nothing;
        name[1] = "Ptent";
        number[1] = Reuse_Ptent;
        name[2] = "full";
        number[2] = Reuse_full;

        setStringToIntegralParameter<int>("MueLu_REUSE", "nothing",
            "Reuse strategy in MueLu contact preconditioner 2.\n"
            "The options are 'none', 'Ptent' and 'full'.\n"
            "'full' means: reuse the full multigrid hierarchy.\n"
            "'Ptent': reuse aggregates and nonsmoothed transfer operator.\n"
            "The MueLu_Reuse parameter only makes sense together with AZREUSE.\n",
            name, number, &list);
      }

      //  // parameters for AMG(BS)
      //  setNumericStringParameter("AMGBS_BS_DAMPING","1.3 1.3 1.3",
      //                            "Relaxation factor for Braess-Sarazin smoother within AMGBS
      //                            method", &list);
      //
      //  setNumericStringParameter("AMGBS_BS_PCSWEEPS","2 2 2",
      //                            "number of jacobi/sgs sweeps for smoothing/solving pressure
      //                            correction equation within Braess-Sarazin. only necessary for
      //                            jacobi/gauss seidel", &list);
      //
      //  setNumericStringParameter("AMGBS_BS_PCDAMPING","1.0 1.0 1.0",
      //                              "jacobi damping factors for smoothing/solving pressure
      //                              correction equation within Braess-Sarazin. only necessary for
      //                              jacobi/gauss seidel", &list);
      //
      //
      //  setStringToIntegralParameter<int>(
      //    "AMGBS_PSMOOTHER_VEL","PA-AMG","Prolongation/Restriction smoothing strategy (velocity
      //    part in AMGBS preconditioner)",
      //    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
      //    tuple<int>(PA_AMG,SA_AMG,PG_AMG,PG2_AMG),
      //    &list);
      //  setStringToIntegralParameter<int>(
      //    "AMGBS_PSMOOTHER_PRE","PA-AMG","Prolongation/Restriction smoothing strategy (pressure
      //    part in AMGBS preconditioner)",
      //    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
      //    tuple<int>(PA_AMG,SA_AMG,PG_AMG,PG2_AMG),
      //    &list);
      //
      //  setStringToIntegralParameter<int>(
      //    "AMGBS_BS_PCCOARSE","Umfpack","approximation algorithm for solving pressure correction
      //    equation (coarsest level)",
      //    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric
      //    Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel
      //    stand-alone"), tuple<int>(0,1,2,3,4,5,6,7,8), &list);
      //
      //  setStringToIntegralParameter<int>(
      //    "AMGBS_BS_PCMEDIUM","Umfpack","approximation algorithm for solving pressure correction
      //    equation (medium level)",
      //    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric
      //    Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel
      //    stand-alone"), tuple<int>(0,1,2,3,4,5,6,7,8), &list);
      //
      //  setStringToIntegralParameter<int>(
      //    "AMGBS_BS_PCFINE","Umfpack","approximation algorithm for solving pressure correction
      //    equation (finest level)",
      //    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric
      //    Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel
      //    stand-alone"), tuple<int>(0,1,2,3,4,5,6,7,8), &list);

      // switch order of blocks in BGS2x2 preconditioner
      setStringToIntegralParameter<int>("BGS2X2_FLIPORDER", "block0_block1_order",
          "BGS2x2 flip order parameter",
          tuple<std::string>("block0_block1_order", "block1_block0_order"), tuple<int>(0, 1),
          &list);

      // damping parameter for BGS2X2
      DoubleParameter(
          "BGS2X2_GLOBAL_DAMPING", 1., "damping parameter for BGS2X2 preconditioner", &list);
      DoubleParameter(
          "BGS2X2_BLOCK1_DAMPING", 1., "damping parameter for BGS2X2 preconditioner block1", &list);
      DoubleParameter(
          "BGS2X2_BLOCK2_DAMPING", 1., "damping parameter for BGS2X2 preconditioner block2", &list);

      // parameters for permutation of linear systems
      {
        Teuchos::Tuple<std::string, 5> name;
        Teuchos::Tuple<int, 5> number;
        name[0] = "none";
        number[0] = Permutation_none;
        name[1] = "Algebraic";
        number[1] = Permutation_algebraic;
        name[2] = "algebraic";
        number[2] = Permutation_algebraic;
        name[3] = "Local";
        number[3] = Permutation_local;
        name[4] = "local";
        number[4] = Permutation_local;

        setStringToIntegralParameter<int>("PERMUTE_SYSTEM", "none",
            "allow linear solver to permute linear system to improve properties of linear system "
            "for iterative methods.",
            name, number, &list);
      }
      DoubleParameter("NON_DIAGDOMINANCE_RATIO", 1.,
          "matrix rows with diagEntry/maxEntry<nonDiagDominanceRatio are marked to be "
          "significantly non-diagonal dominant (default: 1.0 = mark all non-diagonal dominant "
          "rows)",
          &list);

      // verbosity flag (for Belos)
      IntParameter(
          "VERBOSITY", 0, "verbosity level (0=no output,... 10=extreme), for Belos only", &list);

      // the only one stratimikos specific parameter
      StringParameter("STRATIMIKOS_XMLFILE", "none", "xml file for stratimikos parameters", &list);

      // user-given name of solver block (just for beauty)
      StringParameter("NAME", "No_name", "User specified name for solver block", &list);

      // damping parameter for SIMPLE
      DoubleParameter("SIMPLE_DAMPING", 1., "damping parameter for SIMPLE preconditioner", &list);



      // Parameters for AMGnxn Preconditioner
      StringParameter("AMGNXN_TYPE", "AMG(BGS)",
          "Name of the pre-built preconditioner to be used. If set to\"XML\" the preconditioner is "
          "defined using a xml file",
          &list);
      StringParameter(
          "AMGNXN_XML_FILE", "none", "xml file defining the AMGnxn preconditioner", &list);
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

  }  // end of namespace SOLVER
}  // end of namespace INPAR
