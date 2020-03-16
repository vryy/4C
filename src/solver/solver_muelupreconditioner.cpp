/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1

\maintainer Martin Kronbichler
*/
#ifdef HAVE_MueLu

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include <MueLu_ConfigDefs.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>

// MueLu
#include <MueLu.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>

#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>

#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include <MueLu_AggregationExportFactory.hpp>


#include <MueLu_EpetraOperator.hpp>  // Aztec interface

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)
#include "muelu/MueLu_BaciFactoryFactory_decl.hpp"  // Baci specific MueLu factories with xml interface
#endif

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include "solver_muelupreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuPreconditioner::MueLuPreconditioner(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuPreconditioner::replaceAll(
    std::string& str, const std::string& from, const std::string& to)
{
  if (from.empty()) return;
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos)
  {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();  // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    // check whether A is a Epetra_CrsMatrix
    // i.e. no block matrix
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);

    // free old matrix first
    P_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // see whether we use standard ml or our own mlapi operator
    // const bool domuelupreconditioner = mllist_.get<bool>("LINALG::MueLu_Preconditioner",false);

    // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
        Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
        Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

    // remove unsupported flags
    mllist_.remove("aggregation: threshold", false);  // no support for aggregation: threshold TODO

#if 0  // helper routine to export null space vectors (adapt file name by hand)
    std::cout << "*** EXPORT nullspace: BEGIN ***" << std::endl;
    // DO NOT COMMIT THIS STUFF
    // write out nullspace
    std::ofstream os;
    os.open("nspvector.out",std::fstream::trunc);

    Teuchos::RCP<std::vector<double> > nsp = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace");
    int nsdim = mllist_.get<int>("null space: dimension");
    int nEq   = mllist_.get<int>("PDE equations");

    // loop over all nullspace vectors
    for(int nsv = 0; nsv<nsdim; nsv++) {
        os << "VECTORS nsp" << nsv << " float" << std::endl;
        // loop over all nodes
        for (int row = 0; row < nsp->size()/(nsdim); row+=nEq) {
            for (int col = 0; col < nEq; col++) {
                os << std::setprecision(16);
                os << (*nsp)[nsv*nsp->size()/nsdim+row+col];
                os << " ";
            }
            os << std::endl;
        }
    }

    os << std::flush;
    os.close();

    std::cout << "*** EXPORT nullspace: END ***" << std::endl;
#endif

#if 0  // helper routine to export matrix and RHS
    // DO NOT COMMIT THIS STUFF
    LINALG::PrintMatrixInMatlabFormat("F.out",*A,true);

    std::ofstream os2;

    // open file for writing
    os2.open("bF.out",std::fstream::trunc);
    os2 << "%%MatrixMarket matrix array real general" << std::endl;
    os2 << b->Map().NumGlobalElements() << " " << 1 << std::endl;

        int NumMyElements1 = b->Map().NumMyElements();
        int MaxElementSize1 = b->Map().MaxElementSize();
        int* MyGlobalElements1 = b->Map().MyGlobalElements();
        int* FirstPointInElementList1(NULL);
        if (MaxElementSize1!=1) FirstPointInElementList1 = b->Map().FirstPointInElementList();
        double ** A_Pointers = b->Pointers();

        for (int i=0; i<NumMyElements1; i++)
        {
            os2 << std::setw(30) << std::setprecision(16) <<  A_Pointers[0][i];    // print out values of 1. vector (only Epetra_Vector supported, no Multi_Vector)
            os2 << endl;
        }
        os2 << flush;

      // close file
      os2.close();

      dserror("exit");
    // END DO NOT COMMIT THIS STUFF
#endif

    // append user-given factories (for export of aggregates, debug info etc...)
    // Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
    // > aggExpFact = Teuchos::rcp(new
    // MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("aggs_level%LEVELID_proc%PROCID.out",/*UCAggFact.get()*/
    // NULL, /*dropFact.get()*/ NULL)); std::vector<Teuchos::RCP<FactoryBase> > vec;
    // vec.push_back(aggExpFact);


    // read MueLu parameters from XML file (if provided as file name in the STRATIMIKOS_XMLFILE
    // parameter in the solver block of the dat file)
    if (mllist_.isParameter("xml file") && mllist_.get<std::string>("xml file") != "none")
    {
      // use parameters from user-provided XML file

      // xxd -i ${INPUT_FILE_NAME} | sed s/}\;/,0x00}\;/ > ${INPUT_FILE_PATH}.h
      /*#include "test.xml.h"
      std::string test_xml_string (reinterpret_cast<char*>(test_xml));
      std::cout << test_xml_string << std::endl << std::endl;*/



      // use xml file for generating hierarchy
      std::string xmlFileName = mllist_.get<std::string>("xml file");

      // screen output
      if (matrix->Comm().MyPID() == 0)
      {
        Teuchos::RCP<Teuchos::FancyOStream> fos =
            Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
        *fos << "Use XML file " << xmlFileName << " for generating MueLu multigrid hierarchy"
             << std::endl;
      }

      // prepare nullspace vector for MueLu
      int numdf = mllist_.get<int>("PDE equations", -1);
      int dimns = mllist_.get<int>("null space: dimension", -1);
      if (dimns == -1 || numdf == -1)
        dserror("Error: PDE equations or null space dimension wrong.");
      Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector =
          Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
              rowMap, dimns, true);
      Teuchos::RCP<std::vector<double>> nsdata =
          mllist_.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

      for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
      {
        Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
        const size_t myLength = nspVector->getLocalLength();
        for (size_t j = 0; j < myLength; j++)
        {
          nspVectori[j] = (*nsdata)[i * myLength + j];
        }
      }

      mueluOp->SetFixedBlockSize(numdf);

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)
      Teuchos::RCP<MueLu::BaciFactoryFactory<Scalar, GlobalOrdinal, LocalOrdinal, Node>>
          myFactFact = Teuchos::rcp(
              new MueLu::BaciFactoryFactory<Scalar, GlobalOrdinal, LocalOrdinal, Node>());
      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(mueluOp->getRowMap()->getComm()), myFactFact);
#else
      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(mueluOp->getRowMap()->getComm()));
#endif

      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->SetDefaultVerbLevel(MueLu::Extreme);
      H->GetLevel(0)->Set("A", mueluOp);
      H->GetLevel(0)->Set("Nullspace", nspVector);
      H->GetLevel(0)->setlib(Xpetra::UseEpetra);
      H->setlib(Xpetra::UseEpetra);

      // provide extra data from outside
      registerEpetraMapinMueLuLevel(mllist_, std::string("contact masterDofMap"), *H);
      registerEpetraMapinMueLuLevel(mllist_, std::string("contact slaveDofMap"), *H);
      registerEpetraMapinMueLuLevel(mllist_, std::string("contact activeDofMap"), *H);
      registerEpetraMapinMueLuLevel(mllist_, std::string("contact innerDofMap"), *H);
      registerMapinMueLuLevel(mllist_, std::string("non diagonal-dominant row map"), *H);
      registerMapinMueLuLevel(mllist_, std::string("near-zero diagonal row map"), *H);

      mueLuFactory.SetupHierarchy(*H);

      /*Teuchos::RCP<Teuchos::FancyOStream> out = fancyOStream(Teuchos::rcpFromRef(std::cout));
      H->GetLevel(0)->print(*out,MueLu::Extreme);
      H->GetLevel(1)->print(*out,MueLu::Extreme);
      H->GetLevel(2)->print(*out,MueLu::Extreme);*/

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
    else
    {
      // Standard case: use ML parameters from dat file

      // Setup MueLu Hierarchy
      MueLu::MLParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(mllist_ /*, vec*/);
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", mueluOp);
      H->GetLevel(0)->setlib(Xpetra::UseEpetra);
      H->setlib(Xpetra::UseEpetra);
      mueLuFactory.SetupHierarchy(*H);

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }  // if (xml file)
  }    // if (create)
}

#endif
