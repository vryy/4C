/*----------------------------------------------------------------------*/
/*!
\file solver_muelupreconditioner.cpp

\brief Interface class for MueLu preconditioner

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/

/*
 * solver_muelupreconditioner.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: wiesner
 */

#ifdef HAVE_MueLu

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"

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
#ifdef HAVE_Trilinos_Q1_2013
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#endif

#include <MueLu_AggregationExportFactory.hpp>
//#include "muelu_ContactInfoFactory_decl.hpp"


// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

#include "solver_muelupreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuPreconditioner::MueLuPreconditioner( FILE * outfile, Teuchos::ParameterList & mllist )
  : PreconditionerType( outfile ),
    mllist_( mllist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuPreconditioner::Setup( bool create,
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

    // see whether we use standard ml or our own mlapi operator
    //const bool domuelupreconditioner = mllist_.get<bool>("LINALG::MueLu_Preconditioner",false);

    // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC,LO,GO,NO,LMO > > mueluA  = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>(mueluA));

    // remove unsupported flags
    mllist_.remove("aggregation: threshold",false); // no support for aggregation: threshold TODO

#if 0

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

    os << flush;
    os.close();

    std::cout << "*** EXPORT nullspace: END ***" << std::endl;
#endif

#if 0
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
    //Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("aggs_level%LEVELID_proc%PROCID.out",/*UCAggFact.get()*/ NULL, /*dropFact.get()*/ NULL));
    //std::vector<Teuchos::RCP<FactoryBase> > vec;
    //vec.push_back(aggExpFact);

#ifdef HAVE_Trilinos_Q1_2013
    if(mllist_.isParameter("xml file") && mllist_.get<std::string>("xml file") != "none"){
      // use parameters from user-provided XML file

      // use xml file for generating hierarchy
      std::string xmlFileName = mllist_.get<std::string>("xml file");
      std::cout << "Use XML file " << xmlFileName << " for generating MueLu multigrid hierarchy" << std::endl;


      // prepare nullspace vector for MueLu
      int numdf = mllist_.get<int>("PDE equations",-1);
      int dimns = mllist_.get<int>("null space: dimension",-1);
      if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");
      Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = mueluA->getRowMap();

      Teuchos::RCP<MultiVector> nspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
      Teuchos::RCP<std::vector<double> > nsdata = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

      for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
          Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
          const size_t myLength = nspVector->getLocalLength();
          for(size_t j=0; j<myLength; j++) {
                  nspVectori[j] = (*nsdata)[i*myLength+j];
          }
      }

      mueluOp->SetFixedBlockSize(numdf);

      ParameterListInterpreter mueLuFactory(xmlFileName,*(mueluOp->getRowMap()->getComm()));

      Teuchos::RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

      H->SetDefaultVerbLevel(MueLu::Extreme);

      H->GetLevel(0)->Set("A", mueluOp);
      H->GetLevel(0)->Set("Nullspace", nspVector);

      mueLuFactory.SetupHierarchy(*H);

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    } else {
#endif
      // Standard case: use ML parameters from dat file

      // Setup MueLu Hierarchy
      MLParameterListInterpreter mueLuFactory(mllist_/*, vec*/);
      Teuchos::RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", mueluOp);
      mueLuFactory.SetupHierarchy(*H);

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
#ifdef HAVE_Trilinos_Q1_2013
    }
#endif

  }
}

#endif

