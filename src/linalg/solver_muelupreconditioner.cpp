/*
 * solver_muelupreconditioner.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: wiesner
 */

#ifdef HAVE_MueLu

#include "../drt_lib/drt_dserror.H"

/*#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"*/

//#include "linalg_mlapi_operator.H"  // Michael's MLAPI based ML preconditioner
//#include "amgpreconditioner.H"      // Tobias' smoothed aggregation AMG implementation in BACI (only for fluids)

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
//#include <Xpetra_MultiVectorFactory.hpp>

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

#include <MueLu_MLInterpreter_decl.hpp>

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
    const bool domuelupreconditioner = mllist_.get<bool>("LINALG::MueLu_Preconditioner",false);

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC,LO,GO,NO,LMO > > mueluA  = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsOperator<SC,LO,GO,NO,LMO>(mueluA));


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

    // remove unsupported flags
    mllist_.remove("aggregation: threshold",false); // no support for aggregation: threshold TODO

    // Setup MueLu Hierarchy
    Teuchos::RCP<Hierarchy> H = MLInterpreter::Setup(mllist_, mueluOp, nspVector);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    //std::cout << mllist_ << std::endl;

  }
}

#endif

