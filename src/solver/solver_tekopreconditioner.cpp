/*!----------------------------------------------------------------------
\file solver_tekopreconditioner.cpp

\brief Preconditioner class using Trilinos TEKOS framework

\level 2

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

Created on: Jul 13, 2011
*----------------------------------------------------------------------*/

#ifdef HAVE_TEKO

#ifndef HAVE_Stratimikos
#error "Stratimikos is needed by Teko. Please add Stratimikos to your configuration."
#endif

#include "../drt_lib/drt_dserror.H"
#include "solver_tekopreconditioner.H"

#include <Teuchos_TimeMonitor.hpp>

// Teko specific includes
#include <Teko_Utilities.hpp>
#include <Teko_InverseFactory.hpp>
#include <Teko_BlockInvDiagonalStrategy.hpp>
#include <Teko_GaussSeidelPreconditionerFactory.hpp>
#include <Teko_SIMPLEPreconditionerFactory.hpp>

#include "teko/teko_baciepetraoperatorwrapper.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::TekoPreconditioner::TekoPreconditioner( FILE * outfile,
                                                        Teuchos::ParameterList & tekolist)
  : PreconditionerType( outfile ),
    tekolist_( tekolist ),
    type_("undefined")
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

    // put A into a Teuchos::RCP pointer (memory not owned by preconditioner!)
    Pmatrix_ = Teuchos::rcp(A,false);

    // Handles some I/O to the output screen
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    // create BACIEpetraOperatorWrapper
    Teuchos::RCP<LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper> blockA =
        Teuchos::rcp(new LINALG::SOLVER::TEKO::Teko_BACIEpetraOperatorWrapper(Pmatrix_));

    // extract type of preconditioner
    type_ = Params().sublist("Teko Parameters").get<std::string>("Prec Type");

    Teuchos::RCP<Teko::BlockPreconditionerFactory> precFact = Teuchos::null;

    /////////////////////////////////////////////////////// TEKO block Gauss-Seidel
    if(type_=="BGS")
    {
      int numBlocks = blockA->GetBlockRowCount();

      std::vector<Teko::LinearOp > inverseLinearOps;
      inverseLinearOps.reserve(numBlocks);
      for(int i=0; i<numBlocks; i++)
      {
        std::stringstream ssinverse;
        ssinverse << "Inverse" << i+1;
        if(!Params().isSublist(ssinverse.str())) dserror("missing parameter sublists for inverses.");
        if(!Params().sublist(ssinverse.str()).isSublist("Stratimikos Parameters")) dserror("Teko block preconditioners need STRATIMIKOS solver objects for the block inverses.");
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
      precFact = Teuchos::rcp(new Teko::GaussSeidelPreconditionerFactory(Teko::GS_UseLowerTriangle,diagInvStrat));
    }
    /////////////////////////////////////////////////////////////// TEKO SIMPLE
    else if(type_=="SIMPLE")
    {
      if(!Params().isSublist("Inverse1")) dserror("missing parameter sublists for inverse.");
      if(!Params().sublist("Inverse1").isSublist("Stratimikos Parameters")) dserror("Teko block preconditioners need STRATIMIKOS solver objects for the block inverses.");
      Teuchos::ParameterList invParams = Params().sublist("Inverse1").sublist("Stratimikos Parameters");
      std::string invType = invParams.get<std::string>("Linear Solver Type");
      double alpha = Params().sublist("Teko Parameters").get<double>("alpha",0.9);
      Teuchos::RCP<Teko::InverseFactory> inverseFact = Teko::invFactoryFromParamList(invParams,invType);

      // check for second inverse (if available)
      if(Params().isSublist("Inverse2"))
      {
        if(!Params().sublist("Inverse2").isSublist("Stratimikos Parameters")) dserror("Teko block preconditioners need STRATIMIKOS solver objects for the block inverses.");
        Teuchos::ParameterList invParams2 = Params().sublist("Inverse2").sublist("Stratimikos Parameters");
        std::string invType2 = invParams2.get<std::string>("Linear Solver Type");
        Teuchos::RCP<Teko::InverseFactory> inverseFact2 = Teko::invFactoryFromParamList(invParams2,invType2);
        precFact = Teuchos::rcp(new Teko::NS::SIMPLEPreconditionerFactory(inverseFact,inverseFact2,alpha));
      }
      else
        precFact = Teuchos::rcp(new Teko::NS::SIMPLEPreconditionerFactory(inverseFact,alpha));
    }

    if (precFact==Teuchos::null) dserror("unknown type for Teko preconditioner");

    prec_ = Teuchos::rcp(new Teko::Epetra::EpetraBlockPreconditioner(precFact));
    prec_->buildPreconditioner(blockA); // use our Teko_2x2EpetraOperatorWrapper as input
  }
}

#endif /* HAVE_TEKO */
