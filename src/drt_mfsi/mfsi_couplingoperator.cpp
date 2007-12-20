
#ifdef CCADISCRET

#include "mfsi_couplingoperator.H"

#include <Thyra_EpetraThyraWrappers.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::CouplingOperator::CouplingOperator()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::CouplingOperator::CouplingOperator(const Teuchos::RCP<Epetra_Operator> &op,
                                         Teuchos::RCP<const Coupling> domaincoup,
                                         Teuchos::RCP<const Coupling> rangecoup)
{
  initialize(op,domaincoup,rangecoup);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::CouplingOperator::initialize(const Teuchos::RCP<Epetra_Operator> &op,
                                        Teuchos::RCP<const Coupling> domaincoup,
                                        Teuchos::RCP<const Coupling> rangecoup)
{
  op_ = Thyra::nonconstEpetraLinearOp(op);
  domaincoup_ = domaincoup;
  rangecoup_ = rangecoup;

  // We always show the master dof map to the outside world. This is
  // the special case needed for FSI. The linear operator actually
  // expects (or creates) vectors in the slave dof map. (This is why
  // we are here in the first place.) So these are going to be
  // translated.

  if (domaincoup_!=Teuchos::null)
  {
    domainspace_ =
      Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<Scalar> >(
        Thyra::create_VectorSpace(domaincoup_->MasterDofMap()));
  }

  if (rangecoup_!=Teuchos::null)
  {
    rangespace_ =
      Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<Scalar> >(
        Thyra::create_VectorSpace(rangecoup_->MasterDofMap()));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::CouplingOperator::uninitialize(Teuchos::RCP<Epetra_Operator> *op)
{
  if (op) *op = op_->epetra_op();

  op_ = Teuchos::null;
  domaincoup_ = Teuchos::null;
  rangecoup_ = Teuchos::null;

  domainspace_ = Teuchos::null;
  rangespace_ = Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::spmdRange() const
{
  return op_->spmdRange();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::spmdDomain() const
{
  return op_->spmdDomain();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MFSI::CouplingOperator::opSupported(Thyra::ETransp M_trans) const
{
  return op_->opSupported(M_trans);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::ScalarProdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::rangeScalarProdVecSpc() const
{
  if (rangecoup_==Teuchos::null)
  {
    return op_->rangeScalarProdVecSpc();
  }
  else
  {
    return rangespace_;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::ScalarProdVectorSpaceBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::domainScalarProdVecSpc() const
{
  if (domaincoup_==Teuchos::null)
  {
    return op_->domainScalarProdVecSpc();
  }
  else
  {
    return domainspace_;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::CouplingOperator::euclideanApply(
  const Thyra::ETransp M_trans,
  const Thyra::MultiVectorBase<Scalar> &X_in,
  Thyra::MultiVectorBase<Scalar> *Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  if (M_trans!=Thyra::NOTRANS)
    dserror("no transpose supported");

  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > X = Teuchos::rcp(&X_in,false);
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > Y = Teuchos::rcp(Y_inout,false);

  // convert domain vector if needed

  if (domaincoup_!=Teuchos::null)
  {
    Teuchos::RCP<const Epetra_MultiVector> e =
      Thyra::get_Epetra_MultiVector(*domaincoup_->MasterDofMap(), X_in);
    e = domaincoup_->MasterToSlave(e);
    X = Thyra::create_MultiVector(e,domainspace_);
  }

  // prepare to convert range vector if needed

  Teuchos::RCP<Epetra_MultiVector> Yslave;

  if (rangecoup_!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Map> m = rangecoup_->SlaveDofMap();
    Yslave = Teuchos::rcp(new Epetra_Vector(*m));
    //Y = Thyra::create_MultiVector(Yslave,Thyra::create_VectorSpace(m));
    Y = Thyra::create_MultiVector(Yslave,op_->rangeScalarProdVecSpc());
  }

  // actual application

  op_->euclideanApply(M_trans,*X,&*Y,alpha,beta);

  // convert range vector

  if (rangecoup_!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_MultiVector> Ymaster = Thyra::get_Epetra_MultiVector(*rangecoup_->MasterDofMap(), *Y_inout);
    rangecoup_->SlaveToMaster(Yslave,Ymaster);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::LinearOpBase<MFSI::CouplingOperator::Scalar> >
MFSI::CouplingOperator::clone() const
{
  assert(0); // ToDo: Implement when needed
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string MFSI::CouplingOperator::description() const
{
  return op_->description();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::CouplingOperator::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  op_->describe(out,verbLevel);
}


#endif
