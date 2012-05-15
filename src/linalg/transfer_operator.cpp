/*
 * transfer_operator.cpp
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */


#include "transfer_operator.H"
#include "transfer_operator_tentative.H"
#include "transfer_operator_saamg.H"
#include "transfer_operator_pgamg.H"
#include "transfer_operator_pgamg2.H"

LINALG::TransferOperator::TransferOperator(const RCP<SparseMatrix>& A, FILE* outfile) :
 A_(A),
 prolongator_(Teuchos::null),
 restrictor_(Teuchos::null),
 outfile_(outfile)
{

}

// empty destructor
LINALG::TransferOperator::~TransferOperator()
{

}


ostream& LINALG::TransferOperator::Print(std::ostream& os) const
{
  os << Label() << endl;
  return(os);
}

RCP<LINALG::TransferOperator> LINALG::TransferOperatorFactory::Create(const string TransferOperatorType, const Teuchos::RCP<SparseMatrix>& A, FILE* outfile)
{
  try
  {
    LINALG::TransferOperator* op = NULL;

    if(TransferOperatorType == "PA-AMG")
      op = new LINALG::TentativeTransferOperator(A,outfile);
    else if(TransferOperatorType == "SA-AMG")
      op = new LINALG::SAAMGTransferOperator(A,outfile);
    else if(TransferOperatorType == "PG-AMG")
      op = new LINALG::PGAMGTransferOperator(A,outfile);
    else if(TransferOperatorType == "PG2-AMG")
      op = new LINALG::PGAMG2TransferOperator(A,outfile);
    else
        dserror("error");

    return rcp(op);
  }
  catch(string str)
  {
    cout << "Error: TransferOperatorFactory::Create: " << str << endl;
    dserror("upps");
  }
  return null;
}

