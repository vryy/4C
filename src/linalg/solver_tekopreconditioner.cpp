/*
 * solver_tekopreconditioner.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: wiesner
 */

#ifdef TRILINOS_DEV

#include "../drt_lib/drt_dserror.H"
#include "solver_tekopreconditioner.H"

#ifdef TRILINOS_DEV
// Teko specific includes
#include <Teko_Utilities.hpp>
#include <Teko_InverseFactory.hpp>
#include <Teko_GaussSeidelPreconditionerFactory.hpp>

#include "teko_baciepetraoperatorwrapper.H"
#endif

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::TekoPreconditioner::TekoPreconditioner( FILE * outfile,
                                                        Teuchos::ParameterList & tekolist)
  : PreconditionerType( outfile ),
    tekolist_( tekolist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::TekoPreconditioner::Setup( bool create,
                                                Epetra_Operator * matrix,
                                                Epetra_MultiVector * x,
                                                Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    LINALG::BlockSparseMatrixBase* A = dynamic_cast<LINALG::BlockSparseMatrixBase*>( matrix );
    if ( A==NULL )
      dserror( "BlockSparseMatrixBase expected" );

    // free old matrix first
    prec_    = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // put A into a RCP pointer (memory not owned by preconditioner!)
    Pmatrix_ = Teuchos::rcp(A,false);

    // create BACIEpetraOperatorWrapper
    Teuchos::RCP<LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper> blockA =
        Teuchos::rcp(new LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper(Pmatrix_));

    cout << blockA->GetBlockRowCount() << " x " << blockA->GetBlockColCount() << " BlockMatrix" << endl;

    // Handles some I/O to the output screen
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    Teuchos::ParameterList invParams1 = Params().sublist("Primary Inverse").sublist("Stratimikos Parameters");
    std::string invType1 = invParams1.get<std::string>("Linear Solver Type");
    Teuchos::RCP<Teko::InverseFactory> inverseFact1 = Teko::invFactoryFromParamList(invParams1,invType1);
    Teko::InverseLinearOp invOp1 = Teko::buildInverse(*inverseFact1,blockA->GetThyraBlock(0,0));
    invOp1->describe(*out,Teuchos::VERB_DEFAULT);

    Teuchos::ParameterList invParams2 = Params().sublist("Secondary Inverse").sublist("Stratimikos Parameters");
    std::string invType2 = invParams2.get<std::string>("Linear Solver Type");
    Teuchos::RCP<Teko::InverseFactory> inverseFact2 = Teko::invFactoryFromParamList(invParams2,invType2);
    Teko::InverseLinearOp invOp2 = Teko::buildInverse(*inverseFact2,blockA->GetThyraBlock(1,1));
    invOp2->describe(*out,Teuchos::VERB_DEFAULT);

    // build 2x2 block Gauss Seidel preconditioner factory
    Teuchos::RCP<Teko::GaussSeidelPreconditionerFactory> precFact =
        Teuchos::rcp(new Teko::GaussSeidelPreconditionerFactory(Teko::GS_UseLowerTriangle,invOp1,invOp2));
    precFact->describe(*out,Teuchos::VERB_HIGH);
    cout << "Preconditioner factory initialized" << endl;

    prec_ = Teuchos::rcp(new Teko::Epetra::EpetraBlockPreconditioner(precFact));
    cout << "created EpetraBlockPreconditioner" << endl;
    prec_->buildPreconditioner(blockA); // use our Teko_2x2EpetraOperatorWrapper as input
    cout << "Preconditioner built" << endl;


//
//    dserror("check me");
  }
}


#endif /* TRILINOS_DEV */
