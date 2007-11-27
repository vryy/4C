
#ifdef CCADISCRET

#include "../drt_lib/drt_dserror.H"

#include "mfsi_preconditionerfactory.H"

#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_InverseLinearOperator.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::PreconditionerFactory::PreconditionerFactory(
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  structure,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  fluid,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  ale,
  Thyra::ConstLinearOperator<double> const& sfi
  )
  : structure_(structure),
    fluid_(fluid),
    ale_(ale),
    sfidentity_(sfi)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MFSI::PreconditionerFactory::isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc ) const
{
  dserror("MFSI::PreconditionerFactory::isCompatible() not implemented");
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Thyra::PreconditionerBase<double> > MFSI::PreconditionerFactory::createPrec() const
{
  return rcp(new Thyra::DefaultPreconditioner<double>());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::PreconditionerFactory::initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
                                                 Thyra::PreconditionerBase<double> *precOp,
                                                 const Thyra::ESupportSolveUse supportSolveUse) const
{
  Teuchos::RCP< const Thyra::LinearOpBase<double,double> > fsiOp = fwdOpSrc->getOp();
  Teuchos::RCP< const Thyra::DefaultBlockedLinearOp<double> > blockFsiOp =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultBlockedLinearOp<double> >(fsiOp);

  Thyra::ConstLinearOperator<double> structOp = blockFsiOp->getBlock(0,0);
  Thyra::ConstLinearOperator<double> fluidOp = blockFsiOp->getBlock(1,1);
  Thyra::ConstLinearOperator<double> aleOp = blockFsiOp->getBlock(2,2);

  // Build inverse operators
  Thyra::ConstLinearOperator<double> invstruct = inverse(*structure_,structOp,Thyra::IGNORE_SOLVE_FAILURE);
  Thyra::ConstLinearOperator<double> invfluid = inverse(*fluid_,fluidOp,Thyra::IGNORE_SOLVE_FAILURE);
  Thyra::ConstLinearOperator<double> invale = inverse(*ale_,aleOp,Thyra::IGNORE_SOLVE_FAILURE);

  // set special MFSI block preconditioner object to Thyra Wrapper

  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > M =
    Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());
  //Thyra::ConstLinearOperator<double> M = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());
  M->beginBlockFill(5,5);
  M->setBlock(0,0,invstruct.constPtr());
  M->setBlock(1,1,invfluid.constPtr());
  M->setBlock(2,2,invale.constPtr());
  M->setBlock(3,3,sfidentity_.constPtr());
  M->setBlock(4,4,sfidentity_.constPtr());
  M->endBlockFill();

  Teuchos::RCP<Thyra::LinearOpBase<double> > myM = M;
  Thyra::LinearOperator<double> constM = myM;

  Thyra::DefaultPreconditioner<double>* defaultPrec =
    &Teuchos::dyn_cast<Thyra::DefaultPreconditioner<double> >(*precOp);

  (*defaultPrec).initializeRight(constM.constPtr());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::PreconditionerFactory::uninitializePrec(Thyra::PreconditionerBase<double> *prec,
                                                   Teuchos::RCP<const Thyra::LinearOpSourceBase<double> >  *fwdOpSrc,
                                                   Thyra::ESupportSolveUse *supportSolveUse) const
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::PreconditionerFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>&)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> MFSI::PreconditionerFactory::getParameterList()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> MFSI::PreconditionerFactory::unsetParameterList()
{
  return Teuchos::null;
}


#endif
