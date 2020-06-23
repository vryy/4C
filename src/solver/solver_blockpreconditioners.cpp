/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 1

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "../drt_lib/drt_dserror.H"

#include <Xpetra_StridedMap.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_EpetraOperator.hpp>  // Aztec interface

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include "solver_blockpreconditioners.H"

// include header files for concrete implementation
#include "bgs2x2_operator.H"                   // Lena's BGS implementation
#include "solver_cheapsimplepreconditioner.H"  // Tobias' CheapSIMPLE

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuBlockPreconditioner::MueLuBlockPreconditioner(
    FILE* outfile, Teuchos::ParameterList& params)
    : LINALG::SOLVER::PreconditionerType(outfile), params_(params)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  // some typedefs
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  if (create)
  {
    // free old matrix first
    P_ = Teuchos::null;

    // This preconditioner is just to demonstrate how it works
    // The fluid preconditioner shows how to feed MueLu for a 2x2 block system with strided maps
    // The TSI preconditioner demonstrates how it works for a 2x2 block system with blocked maps

    // temporary hack: distinguish between "old" SIMPLER_Operator (for fluid
    // only) and "new" more general test implementation
    // bool mt = params_.get<bool>("MESHTYING",false);
    // bool co = params_.get<bool>("CONTACT",false);
    // bool cstr = params_.get<bool>("CONSTRAINT",false);
    bool fl = params_.isSublist("SIMPLER") ||
              params_.get<bool>(
                  "FLUID", false);  // params_.get<bool>("FLUIDSIMPLE",false); // SIMPLE for fluids
    bool tsi = params_.get<bool>("TSI", false);
    // bool elch = params_.get<bool>("ELCH",false);

    if (fl /*|| elch*/)  // pure fluid problems
    {
      // adapt nullspace for splitted pure fluid problem
      int nv = 0;  // number of velocity dofs
      int np = 0;  // number of pressure dofs
      // int dimns = 0; // number of nullspace vectors (for velocity and pressure dofs)
      int ndofpernode = 0;  // dofs per node
      // int nlnode;

      // const Epetra_Map& fullmap = matrix->OperatorRangeMap();
      // const int length = fullmap.NumMyElements();

      Teuchos::RCP<BlockSparseMatrixBase> A =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
      if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      // fix null space for ML inverses
      Teuchos::ParameterList& inv1 = params_.sublist("SubSmoother1");
      // Teuchos::ParameterList& inv2 = params_.sublist("SubSmoother2");
      ndofpernode = inv1.sublist("NodalBlockInformation").get<int>("number of dofs per node", 0);
      nv = inv1.sublist("NodalBlockInformation").get<int>("number of momentum dofs", 0);
      np = inv1.sublist("NodalBlockInformation").get<int>("number of constraint dofs", 0);
      // dimns = inv1.sublist("NodalBlockInformation").get<int>("nullspace dimension",0);

      // build fluid null space in MueLu format

      if (ndofpernode == 0 || nv == 0 || np == 0)
        dserror("Error: PDE equations or null space dimension wrong.");

      // define strided maps
      std::vector<size_t> stridingInfo;
      stridingInfo.push_back(nv);
      stridingInfo.push_back(np);

      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));

      // create maps
      Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> epetra_fullrangemap =
#ifdef TRILINOS_Q1_2015
          Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));
#else
          Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(Teuchos::rcpFromRef(A->FullRangeMap())));
#endif
      Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> fullrangemap =
          Xpetra::StridedMapFactory<LO, GO, NO>::Build(
              epetra_fullrangemap, stridingInfo, -1 /* stridedBlock */, 0 /*globalOffset*/);
      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap1 =
          Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA11->getRowMap(), stridingInfo,
              xA11->getRowMap()->getIndexBase(), 0 /* stridedBlockId */, 0 /* globalOffset */));
      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap2 =
          Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA22->getRowMap(), stridingInfo,
              xA22->getRowMap()->getIndexBase(), 1 /* stridedBlockId */, 0 /* globalOffset */));
      // Teuchos::RCP<StridedMap> strMap1 = Teuchos::rcp(new StridedMap(xA11->getRowMap(),
      // stridingInfo, 0 /* stridedBlock */, 0 /*globalOffset*/)); Teuchos::RCP<StridedMap> strMap2
      // = Teuchos::rcp(new StridedMap(xA22->getRowMap(), stridingInfo, 1 /* stridedBlock */, 0
      // /*0*/
      // /*globalOffset*/));

      // build map extractor
      std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> xmaps;
      xmaps.push_back(strMap1);
      xmaps.push_back(strMap2);

#ifdef TRILINOS_Q1_2015
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO>> map_extractor =
          Xpetra::MapExtractorFactory<Scalar, LO, GO>::Build(fullrangemap, xmaps);
#else
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, NO>> map_extractor =
          Xpetra::MapExtractorFactory<Scalar, LO, GO, NO>::Build(fullrangemap, xmaps);
#endif

      // build blocked Xpetra operator
      Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp = Teuchos::rcp(
          new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 10));
#ifdef TRILINOS_Q1_2015
      bOp->setMatrix(0, 0, xA11);
      bOp->setMatrix(0, 1, xA12);
      bOp->setMatrix(1, 0, xA21);
      bOp->setMatrix(1, 1, xA22);
#else
      bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11)));
      bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12)));
      bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21)));
      bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22)));
#endif
      bOp->fillComplete();

      // create velocity null space
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector1 =
          Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
              strMap1, nv, true);
      for (int i = 0; i < ndofpernode - 1; ++i)
      {
        Teuchos::ArrayRCP<Scalar> nsValues = nspVector1->getDataNonConst(i);
        int numBlocks = nsValues.size() / (ndofpernode - 1);
        for (int j = 0; j < numBlocks; ++j)
        {
          nsValues[j * (ndofpernode - 1) + i] = 1.0;
        }
      }

      // create pressure null space
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector2 =
          Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
              strMap2, np, true);
      Teuchos::ArrayRCP<Scalar> nsValues2 = nspVector2->getDataNonConst(0);
      for (int j = 0; j < nsValues2.size(); ++j)
      {
        nsValues2[j] = 1.0;
      }

      // use parameters from user-provided XML file
      if (params_.isParameter("xml file") == false ||
          params_.get<std::string>("xml file") == "none")
        dserror(
            "no xml file defined in dat file. MueLu_sym for blocked operators only works with xml "
            "parameters.");

      // use xml file for generating hierarchy
      std::string xmlFileName = params_.get<std::string>("xml file");
      Teuchos::RCP<Teuchos::FancyOStream> fos =
          Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(0);
      *fos << "Use XML file " << xmlFileName << " for generating MueLu multigrid hierarchy"
           << std::endl;

      //////////////////////////////////////// prepare setup
      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(bOp->getRangeMap()->getComm()));


      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

      Teuchos::RCP<MueLu::Level> Finest = H->GetLevel(0);
      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set(
          "A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
                   bOp));

      Finest->Set("Nullspace1", nspVector1);
      Finest->Set("Nullspace2", nspVector2);

      /////////////////////////////////// BEGIN setup

      mueLuFactory.SetupHierarchy(*H);

      ///////////////////////////////////// END setup

      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
    else if (tsi == true)
    {
      Teuchos::RCP<BlockSparseMatrixBase> A =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
      if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> fullrangemap =
#ifdef TRILINOS_Q1_2015
          Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));
#else
          Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(Teuchos::rcpFromRef(A->FullRangeMap())));
#endif

      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));

      // define strided maps
      int nPDE = params_.sublist("Inverse1").get<int>("PDE equations", 1);
      std::vector<size_t> stridingInfo1;
      std::vector<size_t> stridingInfo2;
      stridingInfo1.push_back(nPDE);
      stridingInfo2.push_back(1);

      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap1 =
          Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA11->getRowMap(), stridingInfo1,
              xA11->getRowMap()->getIndexBase(), -1 /* stridedBlock */, 0 /*globalOffset*/));
      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap2 =
          Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA22->getRowMap(), stridingInfo2,
              xA22->getRowMap()->getIndexBase(), -1 /* stridedBlock */, 0 /*globalOffset*/));

      // build map extractor
      std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> xmaps;
      xmaps.push_back(strMap1);
      xmaps.push_back(strMap2);

#ifdef TRILINOS_Q1_2015
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO>> map_extractor =
          Xpetra::MapExtractorFactory<Scalar, LO, GO>::Build(fullrangemap, xmaps);
#else
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, NO>> map_extractor =
          Xpetra::MapExtractorFactory<Scalar, LO, GO, NO>::Build(fullrangemap, xmaps);
#endif

      // build blocked Xpetra operator
      Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp = Teuchos::rcp(
          new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 10));
#ifdef TRILINOS_Q1_2015
      bOp->setMatrix(0, 0, xA11);
      bOp->setMatrix(0, 1, xA12);
      bOp->setMatrix(1, 0, xA21);
      bOp->setMatrix(1, 1, xA22);
#else
      bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11)));
      bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12)));
      bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21)));
      bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22)));
#endif
      bOp->fillComplete();

      ///////////////////////////////////////////////////////////////////////
      // prepare nullspace vectors
      ///////////////////////////////////////////////////////////////////////

      // extract null space for primary field
      int dimns = params_.sublist("Inverse1").get<int>("null space: dimension", 1);
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector1 =
          Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
              xA11->getRowMap(), dimns, true);
      Teuchos::RCP<std::vector<double>> nsdata =
          params_.sublist("Inverse1")
              .get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

      for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
      {
        Teuchos::ArrayRCP<Scalar> nspVector11i = nspVector1->getDataNonConst(i);
        const size_t myLength = nspVector1->getLocalLength();
        for (size_t j = 0; j < myLength; j++)
        {
          nspVector11i[j] = (*nsdata)[i * myLength + j];
        }
      }

      // extract null space for secondary field
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector2 =
          Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
              strMap2, 1, true);
      Teuchos::ArrayRCP<Scalar> nsValues2 = nspVector2->getDataNonConst(0);
      for (int j = 0; j < nsValues2.size(); ++j)
      {
        nsValues2[j] = 1.0;
      }

      // use parameters from user-provided XML file
      if (params_.isParameter("xml file") == false ||
          params_.get<std::string>("xml file") == "none")
        dserror(
            "no xml file defined in dat file. MueLu_sym for blocked operators only works with xml "
            "parameters.");

      // use xml file for generating hierarchy
      std::string xmlFileName = params_.get<std::string>("xml file");
      Teuchos::RCP<Teuchos::FancyOStream> fos =
          Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(0);
      *fos << "Use XML file " << xmlFileName << " for generating MueLu multigrid hierarchy"
           << std::endl;

      //////////////////////////////////////// prepare setup
      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(bOp->getRangeMap()->getComm()));


      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

      Teuchos::RCP<MueLu::Level> Finest = H->GetLevel(0);
      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set(
          "A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
                   bOp));

      Finest->Set("Nullspace1", nspVector1);
      Finest->Set("Nullspace2", nspVector2);

      /////////////////////////////////// BEGIN setup

      mueLuFactory.SetupHierarchy(*H);

      ///////////////////////////////////// END setup

      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
    else
    {
      dserror("MueLuBlockPreconditioner does not support this problem type.");
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::SimplePreconditioner::SimplePreconditioner(
    FILE* outfile, Teuchos::ParameterList& params)
    : LINALG::SOLVER::PreconditionerType(outfile), params_(params)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::SimplePreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    // SIMPLER does not need copy of preconditioning matrix to live
    // SIMPLER does not use the downwinding installed here, it does
    // its own downwinding inside if desired

    // free old matrix first
    P_ = Teuchos::null;

    // temporary hack: distinguish between "old" SIMPLER_Operator (for fluid
    // only) and "new" more general test implementation
    bool mt = params_.get<bool>("MESHTYING", false);
    bool co = params_.get<bool>("CONTACT", false);
    bool cstr = params_.get<bool>("CONSTRAINT", false);
    bool fl = params_.isSublist("SIMPLER") ||
              params_.get<bool>(
                  "FLUID", false);  // params_.get<bool>("FLUIDSIMPLE",false); // SIMPLE for fluids
    bool elch = params_.get<bool>("ELCH", false);
    bool gen = params_.get<bool>("GENERAL", false);

    if (mt || co || cstr)
    {
      // adapt ML null space for contact/meshtying/constraint problems
      Teuchos::RCP<BlockSparseMatrixBase> A =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
      if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      // fix null space for "Inverse1"
      //      {
      //        const Epetra_Map& oldmap = A->FullRowMap();
      //        const Epetra_Map& newmap = A->Matrix(0,0).EpetraMatrix()->RowMap();
      //        LINALG::Solver::FixMLNullspace("Inverse1",oldmap, newmap,
      //        params_.sublist("CheapSIMPLE Parameters").sublist("Inverse1"));
      //      }

      // adapt null space for constraint equations
      // Teuchos::ParameterList& inv2 = params_.sublist("Inverse2");
      Teuchos::ParameterList& inv2 = params_.sublist("CheapSIMPLE Parameters").sublist("Inverse2");
      if (inv2.isSublist("ML Parameters"))
      {
        // Schur complement system (1 degree per "node") -> standard nullspace
        inv2.sublist("ML Parameters").set("PDE equations", 1);
        inv2.sublist("ML Parameters").set("null space: dimension", 1);
        const int plength = (*A)(1, 1).RowMap().NumMyElements();
        Teuchos::RCP<std::vector<double>> pnewns =
            Teuchos::rcp(new std::vector<double>(plength, 1.0));
        // TODO: std::vector<double> has zero length for particular cases (e.g. no Lagrange
        // multiplier on this processor)
        //      -> Teuchos::RCP for the null space is set to NULL in Fedora 12 -> dserror
        //      -> Teuchos::RCP points to a random memory field in Fedora 8 -> Teuchos::RCP for null
        //      space is not NULL
        // Temporary work around (ehrl, 21.12.11):
        // In the case of plength=0 the std::vector<double> is rescaled (size 0 -> size 1, initial
        // value 0) in order to avoid problems with ML (ML expects an Teuchos::RCP for the null
        // space != NULL)
        if (plength == 0) pnewns->resize(1, 0.0);
        inv2.sublist("ML Parameters").set("null space: vectors", &((*pnewns)[0]));
        inv2.sublist("ML Parameters").remove("nullspace", false);
        inv2.sublist("Michael's secret vault")
            .set<Teuchos::RCP<std::vector<double>>>("pressure nullspace", pnewns);
      }

      // P_ = Teuchos::rcp(new
      // LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A,params_.sublist("Inverse1"),params_.sublist("Inverse2"),outfile_));
      P_ = Teuchos::rcp(new LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A,
          params_.sublist("CheapSIMPLE Parameters").sublist("Inverse1"),
          params_.sublist("CheapSIMPLE Parameters").sublist("Inverse2"), outfile_));
    }
    else if (fl || elch)  // CheapSIMPLE for pure fluid problems
    {
      // adapt nullspace for splitted pure fluid problem
      int nv = 0;           // number of velocity dofs
      int np = 0;           // number of pressure dofs
      int ndofpernode = 0;  // dofs per node
      int nlnode;

      const Epetra_Map& fullmap = matrix->OperatorRangeMap();
      const int length = fullmap.NumMyElements();

      Teuchos::RCP<BlockSparseMatrixBase> A =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
      if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      // this is a fix for the old SIMPLER sublist
      // if(!params_.isSublist("Inverse1") && params_.isSublist("SIMPLER"))
      // TODO this if clause can probably go away!
      if (!params_.sublist("CheapSIMPLE Parameters").isSublist("Inverse1") &&
          params_.isSublist("SIMPLER"))
      {
        Teuchos::ParameterList& inv1 =
            params_.sublist("CheapSIMPLE Parameters").sublist("Inverse1");
        inv1 = params_;
        inv1.remove("SIMPLER");
        inv1.remove("Inverse1", false);
        Teuchos::ParameterList& inv2 =
            params_.sublist("CheapSIMPLE Parameters").sublist("Inverse2");
        inv2 = params_.sublist("CheapSIMPLE Parameters").sublist("SIMPLER");
        params_.remove("SIMPLER");
        params_.sublist("CheapSIMPLE Parameters").set("Prec Type", "CheapSIMPLE");
        params_.set("FLUID", true);
      }

      // fix null spae for ML inverses
      // Teuchos::ParameterList& inv1 = params_.sublist("Inverse1");
      Teuchos::ParameterList& inv1 = params_.sublist("CheapSIMPLE Parameters").sublist("Inverse1");
      if (inv1.isSublist("ML Parameters"))
      {
        ndofpernode = inv1.sublist("NodalBlockInformation").get<int>("number of dofs per node", 0);
        nv = inv1.sublist("NodalBlockInformation").get<int>("number of momentum dofs", 0);
        np = inv1.sublist("NodalBlockInformation").get<int>("number of constraint dofs", 0);
        if (ndofpernode == 0) dserror("cannot read numdf from NodalBlockInformation");
        if (nv == 0 || np == 0) dserror("nv or np == 0?");
        nlnode = length / ndofpernode;

        inv1.sublist("ML Parameters").set("PDE equations", nv);
        inv1.sublist("ML Parameters").set("null space: dimension", nv);

        const int vlength = A->Matrix(0, 0).RowMap().NumMyElements();
        Teuchos::RCP<std::vector<double>> vnewns =
            Teuchos::rcp(new std::vector<double>(nv * vlength, 0.0));

        for (int i = 0; i < nlnode; ++i)
        {
          (*vnewns)[i * nv] = 1.0;
          (*vnewns)[vlength + i * nv + 1] = 1.0;
          if (nv > 2) (*vnewns)[2 * vlength + i * nv + 2] = 1.0;
        }
        inv1.sublist("ML Parameters").set("null space: vectors", &((*vnewns)[0]));
        inv1.sublist("ML Parameters").remove("nullspace", false);  // necessary??
        inv1.sublist("Michael's secret vault")
            .set<Teuchos::RCP<std::vector<double>>>("velocity nullspace", vnewns);
      }

      // Teuchos::ParameterList& inv2 = params_.sublist("Inverse2");
      Teuchos::ParameterList& inv2 = params_.sublist("CheapSIMPLE Parameters").sublist("Inverse2");
      if (inv2.isSublist("ML Parameters"))
      {
        inv2.sublist("ML Parameters").set("PDE equations", 1);
        inv2.sublist("ML Parameters").set("null space: dimension", 1);
        const int plength = A->Matrix(1, 1).RowMap().NumMyElements();
        Teuchos::RCP<std::vector<double>> pnewns =
            Teuchos::rcp(new std::vector<double>(plength, 1.0));
        inv2.sublist("ML Parameters").set("null space: vectors", &((*pnewns)[0]));
        inv2.sublist("ML Parameters").remove("nullspace", false);  // necessary?
        inv2.sublist("Michael's secret vault")
            .set<Teuchos::RCP<std::vector<double>>>("pressure nullspace", pnewns);
      }

      P_ = Teuchos::rcp(new LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A,
          params_.sublist("CheapSIMPLE Parameters").sublist("Inverse1"),
          params_.sublist("CheapSIMPLE Parameters").sublist("Inverse2"), outfile_));
    }
    // else if(!params_.isSublist("Inverse1") || !params_.isSublist("Inverse2"))
    else if (gen)  // For a general 2x2 block matrix.  This uses MueLu for AMG, not ML.
    {
      // Remark: we are going to ignore everything which is in the params_ > "CheapSIMPLE
      // Parameters" sublist We need only two sublists, params_ > "Inverse1" and params_ >
      // "Inverse2" containing a "MueLu Parameters" sublist. The "MueLu Parameters" sublist should
      // contain the usual stuff: "xml file","PDE equations","null space: dimension" and "nullspace"


      Teuchos::RCP<BlockSparseMatrixBase> A =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
      if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");


      // Check if we provide everything
      if (not params_.isSublist("Inverse1")) dserror("Inverse1 sublist of params_ not found");
      if (not params_.isSublist("Inverse2")) dserror("Inverse2 sublist of params_ not found");
      Teuchos::ParameterList& sublist1 = params_.sublist("Inverse1");
      Teuchos::ParameterList& sublist2 = params_.sublist("Inverse2");
      if (not sublist1.isSublist("MueLu Parameters"))
        dserror("MueLu Parameters sublist of sublist1 not found");
      else
      {
        Teuchos::ParameterList& MueLuList = sublist1.sublist("MueLu Parameters");
        if (not MueLuList.isParameter("PDE equations"))
          dserror("PDE equations not provided for block 1 of 2");
        if (not MueLuList.isParameter("null space: dimension"))
          dserror("null space: dimension not provided for block 1 of 2");
        if (not MueLuList.isParameter("nullspace"))
          dserror("nullspace not provided for block 1 of 2");
        if (MueLuList.get<std::string>("xml file", "none") == "none")
          dserror("xml file not provided for block 1 of 2");
      }
      if (not sublist2.isSublist("MueLu Parameters"))
        dserror("MueLu Parameters sublist of sublist2 not found");
      else
      {
        Teuchos::ParameterList& MueLuList = sublist1.sublist("MueLu Parameters");
        if (not MueLuList.isParameter("PDE equations"))
          dserror("PDE equations not provided for block 2 of 2");
        if (not MueLuList.isParameter("null space: dimension"))
          dserror("null space: dimension not provided for block 2 of 2");
        if (not MueLuList.isParameter("nullspace"))
          dserror("nullspace not provided for block 2 of 2");
        if (MueLuList.get<std::string>("xml file", "none") == "none")
          dserror("xml file not provided for block 2 of 2");
      }

      P_ = Teuchos::rcp(
          new LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A, sublist1, sublist2, outfile_));
    }
    else
    {
      // cout << "************************************************" << endl;
      // cout << "WARNING: SIMPLE for Fluid? expect bugs..." << endl;
      // cout << "************************************************" << endl;
      // Michaels old CheapSIMPLE for Fluid
      // TODO replace me by CheapSIMPLE_BlockPreconditioner

      // P_ = Teuchos::rcp(new LINALG::SOLVER::SIMPLER_Operator(Teuchos::rcp( matrix, false
      // ),params_,params_.sublist("SIMPLER"),outfile_));
      dserror("old SIMPLE not supported any more");
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BGSPreconditioner::BGSPreconditioner(
    FILE* outfile, Teuchos::ParameterList& params, Teuchos::ParameterList& bgslist)
    : LINALG::SOLVER::PreconditionerType(outfile), params_(params), bgslist_(bgslist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BGSPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    P_ = Teuchos::null;

    int numblocks = bgslist_.get<int>("numblocks");

    if (numblocks == 2)  // BGS2x2
    {
      // check whether sublists for individual block solvers are present
      bool haveprec1 = params_.isSublist("Inverse1");
      bool haveprec2 = params_.isSublist("Inverse2");
      if (!haveprec1 or !haveprec2)
        dserror("individual block solvers for BGS2x2 need to be specified");

      int global_iter = bgslist_.get<int>("global_iter");
      double global_omega = bgslist_.get<double>("global_omega");
      int block1_iter = bgslist_.get<int>("block1_iter");
      double block1_omega = bgslist_.get<double>("block1_omega");
      int block2_iter = bgslist_.get<int>("block2_iter");
      double block2_omega = bgslist_.get<double>("block2_omega");
      bool fliporder = bgslist_.get<bool>("fliporder");

      P_ = Teuchos::rcp(new LINALG::BGS2x2_Operator(Teuchos::rcp(matrix, false),
          params_.sublist("Inverse1"), params_.sublist("Inverse2"), global_iter, global_omega,
          block1_iter, block1_omega, block2_iter, block2_omega, fliporder, outfile_));
    }
    else
      dserror(
          "Block Gauss-Seidel BGS2x2 is currently only implemented for a 2x2 system. Use BGSnxn "
          "for a common block Gauss-Seidel implementation (based on Teko package in Trilinos).");
  }
}
