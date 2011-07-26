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
#include <Teko_BlockInvDiagonalStrategy.hpp>
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

    // Handles some I/O to the output screen
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    int numBlocks = blockA->GetBlockRowCount();

    std::vector<Teko::LinearOp > inverseLinearOps;
    inverseLinearOps.reserve(numBlocks);
    for(int i=0; i<numBlocks; i++)
    {
      std::stringstream ssinverse;
      ssinverse << "Inverse" << i+1;
      if(!Params().isSublist(ssinverse.str())) dserror("missing parameter sublists for inverses.");
      Teuchos::ParameterList invParams = Params().sublist(ssinverse.str()).sublist("Stratimikos Parameters");
      std::string invType = invParams.get<std::string>("Linear Solver Type");
      Teuchos::RCP<Teko::InverseFactory> inverseFact = Teko::invFactoryFromParamList(invParams,invType);
      Teko::InverseLinearOp invOp = Teko::buildInverse(*inverseFact,blockA->GetThyraBlock(i,i));
      inverseLinearOps.push_back(invOp);
    }

    const std::vector<Teko::LinearOp > constInverseLinearOps = inverseLinearOps;

    // build Teko::BlockInvDiagStrategy
    Teuchos::RCP<const Teko::StaticInvDiagStrategy> diagInvStrat = Teuchos::rcp(new Teko::StaticInvDiagStrategy(constInverseLinearOps));

    // build 2x2 block Gauss Seidel preconditioner factory
    Teuchos::RCP<Teko::GaussSeidelPreconditionerFactory> precFact =
        Teuchos::rcp(new Teko::GaussSeidelPreconditionerFactory(Teko::GS_UseLowerTriangle,diagInvStrat));
    //precFact->describe(*out,Teuchos::VERB_HIGH);

    prec_ = Teuchos::rcp(new Teko::Epetra::EpetraBlockPreconditioner(precFact));
    prec_->buildPreconditioner(blockA); // use our Teko_2x2EpetraOperatorWrapper as input
  }
}


#endif /* TRILINOS_DEV */
