
#ifdef CCADISCRET

#include "mfsi_couplingoperator.H"


MFSI::CouplingOperator::CouplingOperator()
{
}


MFSI::CouplingOperator::CouplingOperator(const Teuchos::RCP<Epetra_Operator> &op)
{
  initialize(op);
}


void MFSI::CouplingOperator::initialize(const Teuchos::RCP<Epetra_Operator> &op)
{
  op_ = Thyra::nonconstEpetraLinearOp(op);
}


void MFSI::CouplingOperator::uninitialize(Teuchos::RCP<Epetra_Operator> *op)
{
  if (op) *op = op_->epetra_op();

  op_ = Teuchos::null;
}


Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::spmdRange() const
{
  return op_->spmdRange();
}


Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::spmdDomain() const
{
  return op_->spmdDomain();
}


bool MFSI::CouplingOperator::opSupported(Thyra::ETransp M_trans) const
{
  return op_->opSupported(M_trans);
}


Teuchos::RCP<const Thyra::ScalarProdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::rangeScalarProdVecSpc() const
{
  return op_->rangeScalarProdVecSpc();
}


Teuchos::RCP<const Thyra::ScalarProdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::domainScalarProdVecSpc() const
{
  return op_->domainScalarProdVecSpc();
}


void MFSI::CouplingOperator::euclideanApply(
  const Thyra::ETransp M_trans,
  const Thyra::MultiVectorBase<Scalar> &X_in,
  Thyra::MultiVectorBase<Scalar> *Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  op_->euclideanApply(M_trans,X_in,Y_inout,alpha,beta);
}


Teuchos::RCP<const Thyra::LinearOpBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::clone() const
{
  assert(0); // ToDo: Implement when needed
  return Teuchos::null;
}


std::string MFSI::CouplingOperator::description() const
{
  return op_->description();
}


void MFSI::CouplingOperator::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  op_->describe(out,verbLevel);
}


#endif
