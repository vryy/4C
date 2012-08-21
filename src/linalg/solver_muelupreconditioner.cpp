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
#include "linalg_utils.H"

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
#include <MueLu_UCAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>

#include <MueLu_MLParameterListInterpreter_decl.hpp>

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

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC,LO,GO,NO,LMO > > mueluA  = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsOperator<SC,LO,GO,NO,LMO>(mueluA));

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

    std::ofstream os;

    // open file for writing
    os.open("bF.out",std::fstream::trunc);
    os << "%%MatrixMarket matrix array real general" << std::endl;
    os << b->Map().NumGlobalElements() << " " << 1 << std::endl;

        int NumMyElements1 = b->Map().NumMyElements();
        int MaxElementSize1 = b->Map().MaxElementSize();
        int* MyGlobalElements1 = b->Map().MyGlobalElements();
        int* FirstPointInElementList1(NULL);
        if (MaxElementSize1!=1) FirstPointInElementList1 = b->Map().FirstPointInElementList();
        double ** A_Pointers = b->Pointers();

        for (int i=0; i<NumMyElements1; i++)
        {
            os << std::setw(30) << std::setprecision(16) <<  A_Pointers[0][i];    // print out values of 1. vector (only Epetra_Vector supported, no Multi_Vector)
            os << endl;
        }
        os << flush;

      // close file
      os.close();

      dserror("exit");
    // END DO NOT COMMIT THIS STUFF
#endif

    // append user-given factories (for export of aggregates, debug info etc...)
    //Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("aggs_level%LEVELID_proc%PROCID.out",/*UCAggFact.get()*/ NULL, /*dropFact.get()*/ NULL));
    //Teuchos::RCP<MueLu::ContactInfoFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > contactInfoFact = Teuchos::rcp(new MueLu::ContactInfoFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("info_level%LEVELID_proc%PROCID.vtk",Teuchos::null,Teuchos::null/*nspFact*/));
    //std::vector<Teuchos::RCP<FactoryBase> > vec;
    //vec.push_back(aggExpFact);

    // Setup MueLu Hierarchy
    MLParameterListInterpreter mueLuFactory(mllist_/*, vec*/);
    Teuchos::RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", mueluOp);
    mueLuFactory.SetupHierarchy(*H);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    //std::cout << mllist_ << std::endl;

  }
}

#endif

