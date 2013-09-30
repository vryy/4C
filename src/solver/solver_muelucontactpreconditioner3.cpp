/*
 * solver_muelucontactpreconditioner3.cpp
 *
 *  Created on: Dec 6, 2012
 *      Author: tobias
 */


#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "../drt_lib/drt_dserror.H"

#include <MueLu_ConfigDefs.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
//#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

// MueLu
#include <MueLu.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_SmootherPrototype_decl.hpp>

#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_Aggregates.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_MapTransferFactory.hpp>
#include <MueLu_AggregationExportFactory.hpp>
//#include <MueLu_PermutingSmoother.hpp>
#include <MueLu_IfpackSmoother.hpp>

//#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_AdaptiveSaMLParameterListInterpreter.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

#include "muelu/muelu_ContactAFilterFactory_decl.hpp"
#include "muelu/muelu_ContactTransferFactory_decl.hpp"
#include "muelu/muelu_ContactASlaveDofFilterFactory_decl.hpp"
#include "muelu/MueLu_MyTrilinosSmoother_decl.hpp"

#include "solver_muelucontactpreconditioner3.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactPreconditioner3::MueLuContactPreconditioner3( FILE * outfile, Teuchos::ParameterList & mllist )
  : PreconditionerType( outfile ),
    mllist_( mllist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactPreconditioner3::Setup( bool create,
                                              Epetra_Operator * matrix,
                                              Epetra_MultiVector * x,
                                              Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>( matrix );
    if ( A==NULL )
      dserror( "CrsMatrix expected" );

    // free old matrix first
    P_       = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC,LO,GO,NO,LMO > > mueluA  = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>(mueluA));

    MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Write("fineA.mat", *mueluOp);

    // prepare nullspace vector for MueLu
    int numdf = mllist_.get<int>("PDE equations",-1);
    int dimns = mllist_.get<int>("null space: dimension",-1);
    if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = mueluA->getRowMap();
/*
    Teuchos::RCP<MultiVector> nspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
        Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
        const size_t myLength = nspVector->getLocalLength();
        for(size_t j=0; j<myLength; j++) {
                nspVectori[j] = (*nsdata)[i*myLength+j];
        }
    }
*/
    // remove unsupported flags
    mllist_.remove("aggregation: threshold",false); // no support for aggregation: threshold TODO

    //////////////////////////
    // Setup MueLu Hierarchy

    AdaptiveSaMLParameterListInterpreter mueLuFactory(mllist_/*, vec*/);
    Teuchos::RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->setlib(Xpetra::UseEpetra);
    H->GetLevel(0)->Set("A", mueluOp);
    mueLuFactory.SetupHierarchy(*H);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

  }
}

#endif // HAVE_Trilinos_Q1_2013
#endif // HAVE_MueLu

